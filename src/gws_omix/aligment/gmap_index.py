# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import csv
import json
import os
import re

from gws_core import (File, Folder, InputSpec, IntParam, OutputSpec,
                      TaskInputs, TaskOutputs, Utils, task_decorator)
from gws_core.config.config_types import ConfigParams, ConfigSpecs
from gws_core.io.io_spec import InputSpec, OutputSpec
from gws_core.io.io_spec_helper import InputSpecs, OutputSpecs

from ..base_env.omix_env_task import BaseOmixEnvTask
from ..file.fasta_file import FastaFile


@task_decorator("GmapIndex")
class GmapIndex(BaseOmixEnvTask):
    """
    GmapIndex class.

    Represents a task that wraps GMAP index tool. Mandatory to use GMAP tools.

    Configuration options:
        * `threads`: Multi threading options: number of threads to use. [Default =  8]
    """

    input_specs: InputSpecs = {
        'uncompressed_genome_fasta_file': InputSpec(FastaFile, human_name="", short_description="")
    }
    output_specs: OutputSpecs = {
        'gmap_index_file': OutputSpec(Folder, human_name="", short_description="")
    }
    config_specs: ConfigSpecs = {
        "threads": IntParam(default_value=8, min_value=1, short_description="Number of threads [Default =  8] ")
    }

    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        result_file = Folder()
        result_file = Folder(path=self._output_file_path)
        return {"gmap_index_file": result_file}

    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        thread = params["threads"]
        genome_fasta = inputs["uncompressed_genome_fasta_file"]
        genome_fasta_name = os.path.basename(genome_fasta.path)
        db_name = self._get_output_file_path(genome_fasta_name)
        self._output_file_path = os.path.join(self.working_dir, db_name)

        cmd = [
            "gmap_build -t ", thread,
            " -D ", self._output_file_path,
            "-d ", db_name,
            genome_fasta.path, " 2> tmp.gmap_index.log ; mv tmp.gmap_index.log ",
            self._output_file_path
        ]
        return cmd

    def _get_output_file_path(self, fasta_name):
        return fasta_name + ".gmap_index"
