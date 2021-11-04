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
        annot = inputs["gtf_annotation"]   
        cmd = [
            "gffread -w transcripto.tmp.fa -g ",genome_fasta, annot,
            " ; cat transcripto.tmp.fa  | cut -d " " -f 1 > transcripto.tmp.2.fa ; rm transcripto.tmp.fa ; grep \"^>\" ", genome_fasta,
            " | cut -d \" \" -f 1 > decoys.txt ; sed -i.bak -e 's/>//g' decoys.txt ;  cat transcripto.tmp.2.fa", genome_fasta,
            "  > gentrome.fa.gz  ;  ",
            "salmon index -k 31 -t gentrome.fa.gz -d decoys.txt -p ", thread,
            " -i", self._get_output_file_path(params) ,
            "rm gentrome.fa.gz decoys.txt ;"
        ]
        return cmd
    
    def _get_output_file_path(self, params):
        return params["uncompressed_genome_file"] + ".salmon_index"
