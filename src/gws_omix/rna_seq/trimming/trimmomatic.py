# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com


##################### undooooooooooooooooooooone ##############
import os
import re
from gws_core import (
    ConfigParams, Folder, IntParam, StrParam, Task, InputSpecs, OutputSpecs,
    TaskInputs, TaskOutputs, task_decorator, InputSpec, OutputSpec, ConfigSpecs
)
from gws_omix import FastqFolder

from .trimmomatic_env import TrimommaticShellProxyHelper



@task_decorator("Trimmomatic", human_name="Trimmomatic",
                short_description="Eliminate adapters and low-quality reads",hide=True)
class Trimmomatic(Task):
    """
    Trimmomatic is a software tool that trims and filters high-throughput sequencing data.
    It removes adapter sequences, low-quality reads, and other unwanted artifacts.
    Compatible with both single-end and paired-end data.

    Outputs are named so that, for paired-end data, the final files become:
      <sampleName>_trimmed_{Forward_separator}.fastq.gz
      <sampleName>_trimmed_{Reverse_separator}.fastq.gz
    Then HISAT2 or any other aligner can detect them using the same Forward/Reverse separators.
    """

    input_specs: InputSpecs = InputSpecs({
        'fastq_folder': InputSpec(
            FastqFolder,
            human_name="FASTQ Folder",
            short_description="Folder containing raw FASTQ reads"
        )
    })

    output_specs: OutputSpecs = OutputSpecs({
        'output': OutputSpec(
            FastqFolder,
            human_name="Trimmed Reads",
            short_description="Folder containing trimmed FASTQ files"
        )
    })

    config_specs: ConfigSpecs = ConfigSpecs({
        "threads": IntParam(default_value=2, min_value=2, short_description="Number of threads"),
        "min_len": IntParam(default_value=35, min_value=30, short_description="Minimum read length to keep"),
        "sequencing_type": StrParam(
            allowed_values=["Paired-end", "Single-end"],
            short_description="Choose Single-end or Paired-end"
        ),
        "Forward_separator": StrParam(
            allowed_values=["R1", "1", "r1", " "],
            short_description="Forward read identifier (for Paired-end)"
        ),
        "Reverse_separator": StrParam(
            allowed_values=["R2", "2", "r2", " "],
            short_description="Reverse read identifier (for Paired-end)"
        ),
        "5_prime_hard_trimming_read_size": IntParam(
            default_value=0, min_value=0,
            short_description="Bases to remove from 5' end"
        )
    })

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """Trimmomatic logic, naming final outputs consistently for single/paired reads."""

        # === 1) Retrieve parameters
        fastq_folder: FastqFolder = inputs['fastq_folder']
        threads = params["threads"]
        min_len = params["min_len"]
        seq_type = params["sequencing_type"]
        fwd_sep = params["Forward_separator"]
        rev_sep = params["Reverse_separator"]
        headcrop = params["5_prime_hard_trimming_read_size"]

        shell_proxy = TrimommaticShellProxyHelper.create_proxy(self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, "result")
        os.makedirs(result_path, exist_ok=True)

        # e.g. if fwd_sep="1", forward_pattern => "^(.*?)[-_.](?:1)\w*.fastq.gz$"
        forward_pattern = re.compile(rf"^(.*?)[-_.](?:{re.escape(fwd_sep)})\w*\.fastq\.gz$")
        reverse_pattern = re.compile(rf"^(.*?)[-_.](?:{re.escape(rev_sep)})\w*\.fastq\.gz$")

        if seq_type == "Single-end":
            print("[INFO] Processing Single-end data with Trimmomatic...")

            for fname in os.listdir(fastq_folder.path):
                if fname.endswith(".fastq.gz"):
                    # sample_name is everything before .fastq.gz
                    sample_name = fname.replace(".fastq.gz", "")
                    in_file = os.path.join(fastq_folder.path, fname)
                    out_file = os.path.join(
                        result_path, f"{sample_name}_trimmed.fastq.gz")

                    print(
                        f"[INFO] Trimming single-end: {fname} -> {os.path.basename(out_file)}")

                    trim_cmd = (
                        f"trimmomatic SE -threads {threads} -phred33 {in_file} {out_file} "
                        f"ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 "
                        f"LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:{min_len} HEADCROP:{headcrop}"
                    )
                    print("[DEBUG] Command:", trim_cmd)

                    rc = shell_proxy.run(trim_cmd, shell_mode=True)
                    if rc != 0:
                        raise Exception(f"[ERROR] Trimmomatic failed on single-end file {fname}")

        else:
            print("[INFO] Processing Paired-end data with Trimmomatic...")

            # 1) Identify forward & reverse
            pairs_dict = {}
            for fname in os.listdir(fastq_folder.path):
                if fname.endswith(".fastq.gz"):
                    fwd_match = forward_pattern.match(fname)
                    if fwd_match:
                        sample_name = fwd_match.group(1)
                        print(f"[DEBUG] Found forward read => {fname}, sample => {sample_name}")
                        pairs_dict.setdefault(sample_name, {})['FWD'] = fname
                        continue

                    rev_match = reverse_pattern.match(fname)
                    if rev_match:
                        sample_name = rev_match.group(1)
                        print(f"[DEBUG] Found reverse read => {fname}, sample => {sample_name}")
                        pairs_dict.setdefault(sample_name, {})['REV'] = fname
                        continue

            # 2) Process each matched pair
            for sample_name, files in pairs_dict.items():
                if 'FWD' in files and 'REV' in files:
                    R1_in = os.path.join(fastq_folder.path, files['FWD'])
                    R2_in = os.path.join(fastq_folder.path, files['REV'])

                    R1_out = os.path.join(result_path, f"{sample_name}_trimmed_{fwd_sep}.fastq.gz")
                    R2_out = os.path.join(result_path, f"{sample_name}_trimmed_{rev_sep}.fastq.gz")

                    print(f"[INFO] Trimming paired-end for sample '{sample_name}'")
                    print(f"[DEBUG] Forward => {files['FWD']}, Reverse => {files['REV']}")
                    print(f"[DEBUG] Output => {os.path.basename(R1_out)}, {os.path.basename(R2_out)}")

                    trim_cmd = (
                        f"trimmomatic PE -threads {threads} -phred33 {R1_in} {R2_in} "
                        f"{R1_out} {R2_out} "
                        f"ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:true "
                        f"LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:{min_len} HEADCROP:{headcrop}"
                    )
                    print("[DEBUG] Command:", trim_cmd)

                    rc = shell_proxy.run(trim_cmd, shell_mode=True)
                    if rc != 0:
                        raise Exception(f"[ERROR] Trimmomatic failed on paired-end sample '{sample_name}'")

                else:
                    # if partial pair => skip or raise an error
                    print(f"[WARNING] Sample '{sample_name}' missing forward or reverse, skipping pair.")

        # Final: return the folder with newly trimmed files
        return {
            'output': FastqFolder(result_path)
        }
