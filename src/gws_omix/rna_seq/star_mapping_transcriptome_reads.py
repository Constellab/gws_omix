# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com


import os

from gws_core import (File, InputSpec, IntParam, OutputSpec, TaskInputs,
                      TaskOutputs, task_decorator)
from gws_core.config.config_types import ConfigParams, ConfigSpecs
from gws_core.io.io_spec import InputSpec, OutputSpec
from gws_core.io.io_spec_helper import InputSpecs, OutputSpecs
from gws_omix import FastqFolder

from ..base_env.omix_env_task import BaseOmixEnvTask
from ..file.bam_to_quant_file import BAMToQuantFile
from ..file.star_index_folder import StarIndexFolder


@task_decorator("STARMappingTranscriptomeReads")
class STARMappingTranscriptomeReads(BaseOmixEnvTask):
    """
    STARMappingTranscriptomeReads class. Represents a process that wraps STAR program.

    Mapping is done on the whole genome sequence but only reads overlapping gene feature are kept.

    --> THIS BAM OUTPUT FILE !!! IS !!! COMPATIBLE WITH SALMON QUANTIFICATION TOOL (see SalmonQuantStar class) <---

    Configuration options
        * `threads`: Multi threading options: number of threads to use (min=1, max=7). [Default =  4].
        * `memory`: Memory (RAM in Bytes) usage in Bytes. [Default =  8000000000 bytes (=8GB)].

    """

    input_specs: InputSpecs = {
        'fastq_folder': InputSpec(FastqFolder, human_name="", short_description=""),
        'star_genome_index': InputSpec(StarIndexFolder, human_name="", short_description=""),
        'metadata_file': InputSpec(File, human_name="Metadata File", short_description="Metadata File"),
    }

    output_specs: OutputSpecs = {
        'star_bam_file': OutputSpec(BAMToQuantFile, human_name="", short_description=""),
    }

    config_specs: ConfigSpecs = {
        "threads": IntParam(default_value=2, min_value=2, short_description="Number of threads [Default =  2] "),
        "memory": IntParam(default_value=4000000000, min_value=4000000000, short_description="Memory (RAM in Bytes) usage")
    }

    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        result_file = BAMToQuantFile()
        result_file.path = self._get_output_file_path()
        return {"star_bam_file": result_file}

    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        fq_folder = inputs["fastq_folder"]
        genome_index = inputs["star_genome_index"]
        metadata = inputs["metadata_file"]
        thread = params["threads"]
        ram = params["memory"]

        script_file_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [
            " bash ",
            os.path.join(script_file_dir, "./sh/star_transcriptome_mapping.sh"),
            fq_folder.path,
            thread,
            ram,
            metadata.path,
            genome_index.path
        ]

        return cmd

    def _get_output_file_path(self):
        return os.path.join(self.working_dir, "star_transciptome_mapping")
