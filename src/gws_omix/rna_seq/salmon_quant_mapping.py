# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, IntParam, StrParam, ConfigParams, TaskInputs, TaskOutputs, Utils, Folder
from ..base_env.omix_env_task import BaseOmixEnvTask
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
        'salmon_genome_index': (Folder,)       
    }
    output_specs = { 
        'Salmon_quant_file': (File,) 
    }
    config_specs = {
        "threads": IntParam(default_value=2, min_value=1, short_description="Number of threads [Default =  2] ")
     }
   
    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        result_file = File()
        result_file.path = self._output_file_path
        return {"Salmon_quant_file": result_file}
    
    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        thread = params["threads"]
        reads_forward  = inputs["fastq_files_forward"]
        reads_reverse = inputs["fastq_files_reverse"]
        transcripts_index = inputs["salmon_genome_index"]  
        self._output_file_path = self._get_output_file_path(reads_forward)

        cmd = [
            "salmon quant -i ",transcripts_index,
            " -p ",thread,
            "-l A -1", reads_forward,
            "-2", reads_reverse,
            "--validateMappings -o tmp.dir ",
            " ; mv ./tmp.dir/quant.sf ./  ; rm -rf ./tmp.dir ; mv quant.sf ",  self._output_file_path
        ] 

        return cmd

    def _get_output_file_path(self, fastq_file_name):
        return os.path.join(
            self.working_dir, 
            fastq_file_name + ".Salmon_RNAseq_pseudo_mapping.Salmon_counting.csv"
        )
        