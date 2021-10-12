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
from ..file.gtf_file import GTFFile



@task_decorator("STARIndex")
class STARIndex(BaseOmixEnvTask):
    """
    STAR index class. Represents a process that wraps STAR indexing genome process. This is a mandatory step for RNAseq mapping.
    
    Configuration options
        * `threads`: Multi threading options: number of threads to use (min=1, max=7). [Default =  4].
        * `memory`: Memory (RAM in Bytes) usage in Bytes. (min=8000000000 bytes, 8 GB) [Default =  48000000000 bytes (48 000 000 000 Bytes = 48 GB)].
        * `read_size'`: Reads size minus 1 in bp (ex. 150 bp --> 150 - 1 = 149)". [MANDATORY]
    """

    input_specs = {
        'fasta_genome_sequence': (FastaFile,), 
        'gtf_annotation': (GTFFile,)           
    }
    output_specs = {
        'STAR_index': (File,)
    }
    config_specs = {
        "threads": IntParam(default_value=12, min_value=6, short_description="Number of threads [Default =  12] "),
        "memory": IntParam(default_value=48000000000, min_value=8000000000, short_description="Memory (RAM in Bytes) usage" ),
        "read_size": IntParam(short_description="Reads size minus 1 (in bp; ex.: 150pb reads --> 149)")  
    }
   
    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        result_file = File()
        result_file.path = self.working_dir
        return {"STAR_index": result_file}
 
    
    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        thread = params["threads"]
        ram = params["memory"]
        annotation = params["gtf_annotation"]  
        rd_size = params["read_size"]
        fa_genome = params["fasta_genome_sequence"]                        
       
        cmd = [
            bin_file,
            "--runThreadN",thread,
            "--limitGenomeGenerateRAM",ram,
            "--runMode genomeGenerate --genomeDir",self.working_dir, 
            "--genomeFastaFiles", fa_genome,
            "--sjdbGTFfile",annotation,
            "--sjdbOverhang",rd_size
        ]       
        return cmd
