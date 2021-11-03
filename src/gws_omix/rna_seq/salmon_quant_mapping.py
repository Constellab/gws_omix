# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, IntParam, StrParam, ConfigParams, TaskInputs, TaskOutputs, Utils
from ..base.omix_env_task import BaseOmixEnvTask
from ..file.fasta_file import FastaFile
from ..file.fastq_file import FastqFile

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
        'salmon_genome_index': (File,)       
    }
    output_specs = { 
        'Salmon_quant_file': (File,) 
    }
    config_specs = {
        "threads": IntParam(default_value=12, min_value=1, short_description="Number of threads [Default =  12] ")
     }
   
    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        result_file = File()
        result_file.path = self._get_output_file_path(params)
        return {"Salmon_quant_file": result_file}
    
    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        thread = params["threads"]
        reads_forward  = params["fastq_files_forward"]
        reads_reverse = params["fastq_files_reverse"]
        transcripts_index = params["salmon_genome_index"]  

        cmd = [
            "salmon quant -i ",transcripts_index,
            " -p ",thread,
            "-l A -1", reads_forward,
            "-2", reads_reverse,
            "--validateMappings -o tmp.dir ",
            " ; mv ./tmp.dir/quant.sf ./  ; rm -rf ./tmp.dir ; mv quant.sf ",  self._get_output_file_path(params)
        ] 

        return cmd

    def _get_output_file_path(self, params):
        return params["fastq_files_forward"] + ".Salmon_RNAseq_pseudo_mapping.Salmon_counting"