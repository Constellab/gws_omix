# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import File, IntParam, TaskInputs, TaskOutputs, task_decorator
from gws_core.config.config_types import ConfigParams, ConfigSpecs
from gws_core.io.io_spec import InputSpec, OutputSpec
from gws_core.io.io_spec_helper import InputSpecs, OutputSpecs

from ..base_env.omix_env_task import BaseOmixEnvTask
from ..file.fasta_file import FastaFile
from ..file.gtf_file import GTFFile
from ..file.star_index_folder import StarIndexFolder


@task_decorator("STARIndex")
class STARIndex(BaseOmixEnvTask):
    """
    STAR index class. Represents a process that wraps STAR indexing genome process. This is a mandatory step for RNAseq mapping.

    Configuration options
        * `threads`: Multi threading options: number of threads to use (min=1, max=7). [Default =  4].
        * `memory`: Memory (RAM in Bytes) usage in Bytes. (min=8000000000 bytes, 8 GB) [Default =  48000000000 bytes (48 000 000 000 Bytes = 48 GB)].
        * `read_size'`: Reads size minus 1 in bp (ex. 150 bp --> 150 - 1 = 149)". [MANDATORY]
    """

    input_specs: InputSpecs = {
        'genome_file': InputSpec(FastaFile, short_description="Fasta file", human_name="Fasta file"),
        'annotation_file': InputSpec(GTFFile, short_description="GTF file", human_name="GTF annotation file")
    }
    output_specs: OutputSpecs = {
        'star_index': OutputSpec(StarIndexFolder,)
    }
    config_specs: ConfigSpecs = {
        "threads": IntParam(default_value=2, min_value=2, short_description="Number of threads [Default =  2] "),
        "memory": IntParam(default_value=2000000000, min_value=2000000000, short_description="Memory (RAM in Bytes) usage"),
        "read_size": IntParam(short_description="Reads size minus 1 (in bp; ex.: 150pb reads --> 149)")
    }

    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        result_folder = StarIndexFolder()
        result_folder.path = self._get_output_file_path()
        return {"star_index": result_folder}

    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        genome_fasta_file = inputs["genome_file"]
        annotation_gtf_file = inputs["annotation_file"]
        thread = params["threads"]
        ram_mem = params["memory"]
        reads_size = params["read_size"]
        script_file_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [
            " bash ",
            os.path.join(script_file_dir, "./sh/star_index_cmd.sh"),
            thread,
            ram_mem,
            genome_fasta_file.path,
            annotation_gtf_file.path,
            reads_size
        ]

        return cmd

    def _get_output_file_path(self):
        return os.path.join(self.working_dir, "indexed_genome_folder")
