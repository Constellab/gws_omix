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
from ..file.deepec_file import DeepEcFile

@task_decorator("DeepEC")
class DeepEC(DeepECEnvTask):
    """
    DeepEC class. Represents a process that wraps DeepEC program.

    """

    input_specs = {
        'fasta_file': (FastaFile,)
    }
    output_specs = {
        "deepec_output_file": (DeepEcFile,)
    }
    
    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        result_file = File()
        result_file.path = self._get_output_file_path(params)
        return {"deepec_output_file": result_file} 

    
    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:   
        fa_file = params["fasta"]
        cmd = [ 
        " deepec.py -i", fa_file, 
        "-o ", self.working_dir,
        " ; mv DeepEC_Result.txt ", self._get_output_file_path(params)
        ]           
        return cmd

    def _get_output_file_path(self, params):
        return params["fasta_file"]+ ".csv"