# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, IntParam, ConfigParams, TaskInputs, TaskOutputs, Utils, Folder
from ..base_env.omix_env_task import BaseOmixEnvTask
from ..file.fasta_file import FastaFile
from ..file.gtf_file import GTFFile

@task_decorator("SalmonIndex")
class SalmonIndex(BaseOmixEnvTask):
    """
    SalmonIndex class. Represents a process that wraps Salmon index tool. Mandatory to use Salmon pseudo-aligment RNAseq mapping tool (wraps in SalmonQuantMapping class).
    
    Configuration options
        * `threads`: Multi threading options: number of threads to use. [Default =  12]
    """

    input_specs = {
        'uncompressed_genome_file': (FastaFile,),
        'gtf_annotation': (GTFFile,)          
    }
    output_specs = {
        'salmon_index_folder': (Folder,)
    }
    config_specs = {
        "threads": IntParam(default_value=12, min_value=1, short_description="Number of threads [Default =  12] ")
    }
   
    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        result_file = File()
        result_file.path = self._get_output_file_path(params)
        return {"salmon_index_folder": result_file} 
    
    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        thread = params["threads"]
        genome_fasta = inputs["uncompressed_genome_file"]
        genome_fasta_file_name = os.path.basename(genome_fasta.path)
        annot = inputs["gtf_annotation"]   
        self._output_file_path = self._get_output_file_path(genome_fasta_file_name)

        script_file_dir = os.path.dirname(os.path.realpath(__file__))

        cmd = [
            "bash",
            os.path.join(script_file_dir, "./sh/salmon_index_cmd.sh"),
            genome_fasta.path,
            annot.path,
            thread,
            self._output_file_path
        ]
        return cmd
    
    def _get_output_file_path(self, fasta_file_name):
        return os.path.join(
            self.working_dir, 
            fasta_file_name + ".salmon_index"
        )
