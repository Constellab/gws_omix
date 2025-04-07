# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import re
from gws_core import (
    ConfigParams, Folder, IntParam, StrParam, Task, InputSpecs, OutputSpecs,
    TaskInputs, TaskOutputs, task_decorator, InputSpec, OutputSpec, ConfigSpecs
)
from gws_omix import FastqFolder

# We reuse the same ShellProxyHelper class for convenience,
# but it can be renamed or replaced as needed for your environment.
from .fastp_env import FastpShellProxyHelper

@task_decorator("Fastp", human_name="Fastp",
                short_description="Eliminate adapters and low-quality reads using fastp")
class Fastp(Task):
    """
    Fastp is a tool designed to provide fast all-in-one preprocessing for FastQ files.
    It performs adapter trimming, quality filtering, per-read quality pruning (like sliding windows),
    and can be used for both single-end and paired-end reads.

    - Single-end: trims all files ending with .fastq.gz.
    - Paired-end: uses regex to detect forward vs. reverse read files and
                  names outputs using user-supplied separators.
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

    config_specs: ConfigSpecs = {
        "threads": IntParam(default_value=2, min_value=2,
                            short_description="Number of threads"),
        "sequencing_type": StrParam(allowed_values=["Paired-end", "Single-end"],
                                    short_description="Choose Single-end or Paired-end"),
        "Forward_separator": StrParam(allowed_values=["R1", "1", "r1", " "],
                                      short_description="Forward read identifier (for Paired-end)"),
        "Reverse_separator": StrParam(allowed_values=["R2", "2", "r2", " "],
                                      short_description="Reverse read identifier (for Paired-end)"),
        "5_prime_hard_trimming_read_size": IntParam(default_value=0, min_value=0,
                                                    short_description="Bases to remove from 5' end")
    }

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        # === 1) Retrieve user parameters
        fastq_folder: FastqFolder = inputs['fastq_folder']
        threads = params["threads"]
        seq_type = params["sequencing_type"]
        fwd_sep = params["Forward_separator"]
        rev_sep = params["Reverse_separator"]
        headcrop = params["5_prime_hard_trimming_read_size"]

        # === 2) Create working directory
        shell_proxy = FastpShellProxyHelper.create_proxy(self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, "result")
        os.makedirs(result_path, exist_ok=True)

        # === 3) Build dynamic regex for forward/reverse
        # Example:
        #   if fwd_sep="r1", forward_pattern="^(.*?)[-_.](?:r1)\w*.fastq.gz$"
        #   group(1) = sample name
        forward_pattern = re.compile(rf"^(.*?)[-_.](?:{re.escape(fwd_sep)})\w*\.fastq\.gz$")
        reverse_pattern = re.compile(rf"^(.*?)[-_.](?:{re.escape(rev_sep)})\w*\.fastq\.gz$")

        # === 4) Single-end logic
        if seq_type == "Single-end":
            print("[INFO] Processing Single-end data with fastp...")
            for fname in os.listdir(fastq_folder.path):
                if fname.endswith(".fastq.gz"):
                    sample_name = fname.replace(".fastq.gz", "")
                    in_file = os.path.join(fastq_folder.path, fname)
                    out_file = os.path.join(result_path, f"{sample_name}_trimmed.fastq.gz")

                    print(f"[INFO] Trimming single-end: {fname} -> {os.path.basename(out_file)}")

                    # Construct fastp command for single-end
                    fastp_cmd = (
                        f"fastp "
                        f"--in1 {in_file} "
                        f"--out1 {out_file} "
                        f"--trim_front1 {headcrop} "
                        f"--thread {threads} "
                    )
                    print("[DEBUG] Command:", fastp_cmd)

                    rc = shell_proxy.run(fastp_cmd, shell_mode=True)
                    if rc != 0:
                        raise Exception(f"fastp failed on single-end file {fname}")

        # === 5) Paired-end logic
        else:
            print("[INFO] Processing Paired-end data with fastp...")

            # Step 1: Identify forward & reverse reads
            pairs_dict = {}
            for fname in os.listdir(fastq_folder.path):
                if fname.endswith(".fastq.gz"):
                    fwd_match = forward_pattern.match(fname)
                    if fwd_match:
                        sample_name = fwd_match.group(1)
                        pairs_dict.setdefault(sample_name, {})['FWD'] = fname
                        continue

                    rev_match = reverse_pattern.match(fname)
                    if rev_match:
                        sample_name = rev_match.group(1)
                        pairs_dict.setdefault(sample_name, {})['REV'] = fname
                        continue

            # Step 2: Run fastp on each matched pair
            for sample_name, files in pairs_dict.items():
                if 'FWD' in files and 'REV' in files:
                    R1_in = os.path.join(fastq_folder.path, files['FWD'])
                    R2_in = os.path.join(fastq_folder.path, files['REV'])

                    # Output filenames
                    R1_out = os.path.join(result_path, f"{sample_name}_trimmed_{fwd_sep}.fastq.gz")
                    R2_out = os.path.join(result_path, f"{sample_name}_trimmed_{rev_sep}.fastq.gz")

                    print(f"[INFO] Trimming paired reads: {files['FWD']} & {files['REV']}")
                    print(f"[DEBUG] Outputs -> {os.path.basename(R1_out)}, {os.path.basename(R2_out)}")

                    fastp_cmd = (
                        f"fastp "
                        f"--in1 {R1_in} --in2 {R2_in} "
                        f"--out1 {R1_out} --out2 {R2_out} "
                        f"--trim_front1 {headcrop} --trim_front2 {headcrop} "
                        f"--thread {threads} "
                        f"--detect_adapter_for_pe"
                    )
                    print("[DEBUG] Command:", fastp_cmd)

                    rc = shell_proxy.run(fastp_cmd, shell_mode=True)
                    if rc != 0:
                        raise Exception(f"fastp failed on paired-end sample {sample_name}")

        # === 6) Return the result folder
        return {'output': FastqFolder(result_path)}
