# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, IntParam, ConfigParams, TaskInputs, TaskOutputs, Utils
from ..base.omix_env_task import BaseOmixEnvTask
from ..file.fasta_file import FastaFile

@task_decorator("GmapIndex")
class GmapIndex(BaseOmixEnvTask):
    """
    GMAP index class. Represents a process that wraps GMAP index tool. Mandatory to use GMAP tools.
    
    Configuration options
        * `threads`: Multi threading options: number of threads to use. [Default =  8]
    """

    input_specs = {
        'uncompressed_genome_fasta_file': (FastaFile,),      
    }
    output_specs = {
        'Gmap_index_file': (File,)
    }
    config_specs = {
        "threads": IntParam(default_value=8, min_value=1, description="Number of threads [Default =  8] ")
    }
   
    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        result_file = File()
        result_file.path = self._get_output_file_path(params)
        return {"salmon_index_file": result_file} 
    
    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
         
        thread = params["threads"]
        genome_fasta = params["uncompressed_genome_fasta_file"]
 
        cmd = [
            "gmap_build -t ", thread,
            " -D ", self.working_dir,
            "-d ", self._get_output_file_path(params),
            genome_fasta
        ]           

        return cmd

    def _get_output_file_path(self, params):
        return params["uncompressed_genome_fasta_file"] + ".gmap_index"
