# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, IntParam, TaskInputs, TaskOutputs,
                      task_decorator)

from ..base_env.omix_env_task import BaseOmixEnvTask
from ..file.fastq_file import FastqFile
from ..file.salmon_index_result_folder import SalmonIndexResultFolder
from ..file.salmon_quant_mapping_output_file import SalmonQuantMappingFile


@task_decorator("SalmonQuantMapping")
class SalmonQuantMapping(BaseOmixEnvTask):
    """
    SalmonQuantMapping class. Represents a process that wraps Salmon pseudo-aligment RNAseq mapping program (salmon quant).

    Configuration options
        * `threads`: Multi threading options: number of threads to use (min=6). [Default =  12].
    """

    input_specs = {
        'fastq_files_forward': (FastqFile,),
        'fastq_files_reverse': (FastqFile,),
        'salmon_genome_index': (SalmonIndexResultFolder,)
    }
    output_specs = {
        'salmon_quant_file': (SalmonQuantMappingFile,)
    }
    config_specs = {
        "threads": IntParam(default_value=2, min_value=1, short_description="Number of threads [Default =  2] ")
    }

    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        path = self._get_output_file_path(inputs)
        result_file = SalmonQuantMappingFile(path=path)
        return {"salmon_quant_file": result_file}

    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        thread = params["threads"]
        reads_forward = inputs["fastq_files_forward"]
        reads_reverse = inputs["fastq_files_reverse"]
        transcripts_index = inputs["salmon_genome_index"]

        cmd = [
            "salmon quant -i ", transcripts_index.path,
            " -p ", thread,
            "-l A -1", reads_forward.path,
            "-2", reads_reverse.path,
            "--validateMappings -o ./ ; ",
            " mv quant.sf ",  self._get_fastq_file_name(inputs)
        ]

        return cmd

    def _get_fastq_file_name(self, inputs):
        fastq = inputs["fastq_files_forward"]
        return os.path.basename(fastq.path) + ".Salmon_counting.csv"

    def _get_output_file_path(self, inputs):
        return os.path.join(
            self.working_dir,
            self._get_fastq_file_name(inputs)
        )
