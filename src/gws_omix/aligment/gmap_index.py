# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, IntParam, ConfigParams, TaskInputs, TaskOutputs, Utils, Folder
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
        'gmap_index_file': (Folder,)
    }
    config_specs = {
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
        return  fasta_name + ".gmap_index"
