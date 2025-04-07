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
from .hisat2_env import Hisat2ShellProxyHelper


@task_decorator("Hisat2_Align", human_name="Hisat2_Align",
                short_description="Align RNA-seq reads to a reference genome using HISAT2")
class Hisat2Align(Task):
    """
    This task performs RNA-seq read alignment using HISAT2 (a splice-aware aligner).
    It supports both single-end and paired-end reads, and outputs sorted BAM files.
    """

    input_specs: InputSpecs = InputSpecs({
        'fastq_folder': InputSpec(FastqFolder, human_name="Trimmed Fastq folder",
                                  short_description="Folder containing trimmed FASTQ reads"),
        'hisat2_index': InputSpec(Folder, human_name="HISAT2 Index",
                                  short_description="Folder containing HISAT2 genome index")
    })

    output_specs: OutputSpecs = OutputSpecs({
        'output': OutputSpec(Folder, human_name="Aligned BAM files",
                             short_description="Folder containing sorted BAM files")
    })

    config_specs: ConfigSpecs = {
        "threads": IntParam(default_value=4, min_value=2, short_description="Number of threads"),
        "sequencing_type": StrParam(allowed_values=["Paired-end", "Single-end"],
                                    short_description="Sequencing type"),
        "Forward_separator": StrParam(allowed_values=["R1", "1", "r1", " "],
                                      short_description="Forward read identifier for Paired-end"),
        "Reverse_separator": StrParam(allowed_values=["R2", "2", "r2", " "],
                                      short_description="Reverse read identifier for Paired-end")
    }

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """Main entry point: run the HISAT2 alignment."""

        # Retrieve user params
        hisat2_index: Folder = inputs['hisat2_index']
        fastq_folder: FastqFolder = inputs['fastq_folder']
        threads = params["threads"]
        seq_type = params["sequencing_type"]
        fwd_sep = params["Forward_separator"]
        rev_sep = params["Reverse_separator"]

        # Assume index prefix is "genome" inside hisat2_index
        hisat2_index_path = os.path.join(hisat2_index.path, "genome_index")

        # Create output directory
        shell_proxy = Hisat2ShellProxyHelper.create_proxy(self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'result')
        os.makedirs(result_path, exist_ok=True)

        # Regex to detect forward & reverse files, e.g. ".*?[-_.](?:R1|1|r1)\w*\.fastq.gz$"
        forward_pattern = re.compile(rf"^(.*?)[-_.](?:{re.escape(fwd_sep)})\w*\.fastq\.gz$")
        reverse_pattern = re.compile(rf"^(.*?)[-_.](?:{re.escape(rev_sep)})\w*\.fastq\.gz$")

        # == Single-end logic ==
        if seq_type == "Single-end":
            print("[INFO] Single-end alignment with HISAT2...")
            for fname in os.listdir(fastq_folder.path):
                if fname.endswith(".fastq.gz"):
                    sample_name = fname.replace(".fastq.gz", "")
                    fastq_path = os.path.join(fastq_folder.path, fname)

                    # Output sorted BAM
                    bam_output = os.path.join(result_path, f"{sample_name}.bam")

                    hisat2_cmd = (
                        f"hisat2 -q --rna-strandness F -p {threads} "
                        f"-x {hisat2_index_path} -U {fastq_path} | "
                        f"samtools sort -o {bam_output}"
                    )
                    print("[DEBUG] Single-end CMD:", hisat2_cmd)

                    rc = shell_proxy.run(hisat2_cmd, shell_mode=True)
                    if rc != 0:
                        raise Exception(f"Error during HISAT2 alignment (single-end) for {sample_name}")

        # == Paired-end logic ==
        else:
            print("[INFO] Paired-end alignment with HISAT2...")

            # 1) Identify forward & reverse read files
            pairs_dict = {}

            for fname in os.listdir(fastq_folder.path):
                if fname.endswith(".fastq.gz"):
                    # Check forward
                    fwd_m = forward_pattern.match(fname)
                    if fwd_m:
                        sample_name = fwd_m.group(1)
                        pairs_dict.setdefault(sample_name, {})['FWD'] = fname
                        continue

                    # Check reverse
                    rev_m = reverse_pattern.match(fname)
                    if rev_m:
                        sample_name = rev_m.group(1)
                        pairs_dict.setdefault(sample_name, {})['REV'] = fname
                        continue

            # 2) Run HISAT2 for each matched pair
            for sample_name, files in pairs_dict.items():
                if 'FWD' in files and 'REV' in files:
                    R1_in = os.path.join(fastq_folder.path, files['FWD'])
                    R2_in = os.path.join(fastq_folder.path, files['REV'])

                    bam_output = os.path.join(result_path, f"{sample_name}.bam")

                    hisat2_cmd = (
                        f"hisat2 -q --rna-strandness FR -p {threads} "
                        f"-x {hisat2_index_path} -1 {R1_in} -2 {R2_in} | "
                        f"samtools sort -o {bam_output}"
                    )
                    print("[DEBUG] Paired-end CMD:", hisat2_cmd)

                    rc = shell_proxy.run(hisat2_cmd, shell_mode=True)
                    if rc != 0:
                        raise Exception(f"Error during HISAT2 alignment (paired-end) for {sample_name}")

        return {'output': Folder(result_path)}
