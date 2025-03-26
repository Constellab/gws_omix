# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import re
from gws_core import (ConfigParams, Folder, InputSpec, IntParam, OutputSpec, InputSpecs, OutputSpecs,
                      Task, TaskInputs, TaskOutputs, task_decorator, ConfigSpecs, StrParam)

from gws_omix import FastqFolder

from .fastq_init_env import FastqInitShellProxyHelper


@task_decorator("FastQC", human_name="FastQC",
                short_description="Check the quality of reads using fastqc")
class FastqcInit(Task):
    """
    FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines.
    It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.
    """

    input_specs: InputSpecs = InputSpecs({
        'fastq_folder': InputSpec(FastqFolder, human_name="Fastq folder",
                                  short_description="folder containing reads")})
    output_specs: OutputSpecs = OutputSpecs({
        'output': OutputSpec(Folder, human_name="Quality result",
                             short_description="reads quality report")})

    config_specs: ConfigSpecs = ConfigSpecs({
        "threads": IntParam(default_value=8, min_value=2, short_description="Number of threads"),
        "sequencing_mode": StrParam(default_value="paired", allowed_values=["paired", "single"],
                                    short_description="Sequencing mode: paired-end or single-end")
    })

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # Retrieve input parameters
        input_folder: FastqFolder = inputs['fastq_folder']
        thrd = params["threads"]
        sequencing_mode = params["sequencing_mode"]

        # Create working directory
        shell_proxy = FastqInitShellProxyHelper.create_proxy(
            self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'result')
        os.makedirs(result_path)

        # Regular expression for detecting FASTQ files
        r1_file_pattern = re.compile(r'.*?[-_.][R1|1|r1]\w*\.fastq\.gz')
        r2_file_pattern = re.compile(r'.*?[-_.][R2|2|r2]\w*\.fastq\.gz')

        r1_files = [f for f in os.listdir(
            input_folder.path) if r1_file_pattern.match(f)]
        r2_files = [f for f in os.listdir(
            input_folder.path) if r2_file_pattern.match(f)]

        if sequencing_mode == "single":
            # In single-end mode, include all fastq.gz files (even those without R1 or 1)
            all_fastq_files = [f for f in os.listdir(
                input_folder.path) if f.endswith(".fastq.gz")]
            fastqc_cmd = f'fastqc {" ".join([os.path.join(input_folder.path, f) for f in all_fastq_files])} -o {result_path} -t {thrd}'
        else:
            # In paired-end mode, both R1 and R2 files are considered
            fastqc_cmd = f'fastqc {" ".join([os.path.join(input_folder.path, f) for f in r1_files])} {" ".join([os.path.join(input_folder.path, f) for f in r2_files])} -o {result_path} -t {thrd}'

        res = shell_proxy.run(fastqc_cmd, shell_mode=True)

        if res != 0:
            raise Exception("An error occurred when formatting output files")

        folder = Folder(result_path)

        # Return the output
        return {'output': folder}
