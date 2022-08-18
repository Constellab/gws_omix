# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import csv
import json
import os
import re

from gws_core import (File, InputSpec, OutputSpec, TaskInputs, TaskOutputs,
                      Utils, task_decorator)
from gws_core.config.config_types import ConfigParams, ConfigSpecs
from gws_core.io.io_spec import InputSpec, OutputSpec
from gws_core.io.io_spec_helper import InputSpecs, OutputSpecs

from ..base_env.deepec_env_task import DeepECEnvTask
from ..file.deepec_file import DeepECFile
from ..file.fasta_file import FastaFile


@task_decorator("DeepEC")
class DeepEC(DeepECEnvTask):
    """
    DeepEC class. Represents a process that wraps DeepEC program.

    """

    input_specs: InputSpecs = {
        'fasta_file': InputSpec(FastaFile, human_name="", short_description="")
    }
    output_specs: OutputSpecs = {
        'deepec_file': OutputSpec(DeepECFile, human_name="", short_description="")
    }

    def gather_outputs(self,  inputs: TaskInputs) -> TaskOutputs:
        result_file = DeepECFile()
        result_file.path = self._output_file_path
        return {"deepec_file": result_file}

    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        fasta_file = inputs["fasta_file"]
        fasta_file_name = os.path.basename(fasta_file.path)
        self._output_file_path = self._get_output_file_path(fasta_file_name)

        script_file_dir = os.path.dirname(os.path.realpath(__file__))
        deepec_py = "/lab/.sys/lib/deepec/deepec.py"  # TODO use a variable
        cmd = [
            " bash ",
            os.path.join(script_file_dir, "./sh/deepec_cmd.sh"),
            deepec_py,
            fasta_file.path,
            self.working_dir,
            self._output_file_path
        ]

        return cmd

    def _get_output_file_path(self, fasta_file_name):
        return os.path.join(
            self.working_dir,
            fasta_file_name + ".deepec_output"
        )
