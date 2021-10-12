# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, ConfigParams, TaskInputs, TaskOutputs, Utils
from ..base.deepec_env_task import DeepECEnvTask
from ..file.fasta_file import FastaFile
from ..file.deepec_file import DeepECFile

@task_decorator("DeepEC")
class DeepEC(DeepECEnvTask):
    """
    DeepEC class. Represents a process that wraps DeepEC program.

    """

    input_specs = {
        'fasta_file': (FastaFile,)
    }
    output_specs = {
        'deepec_file': (DeepECFile,)
    }
    
    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        result_file = DeepECFile()
        result_file.path = self._output_file_path
        return {"deepec_file": result_file} 
    
    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:   
        fasta_file = inputs["fasta_file"]
        fasta_file_name = os.path.basename(fasta_file.path)
        self._output_file_path = self._get_output_file_path(fasta_file_name)

        script_file_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [ 
            " bash ", 
            os.path.join(script_file_dir, "./sh/deepec_cmd.sh"),         
            fasta_file.path, 
            self.working_dir,
            self._output_file_path
        ]
          
        return cmd

    def _get_output_file_path(self, fasta_file_name) :
        return os.path.join(
            self.working_dir, 
            fasta_file_name + ".deepec_output"
        )
