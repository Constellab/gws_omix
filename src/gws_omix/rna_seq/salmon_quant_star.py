# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, IntParam, StrParam, ConfigParams, TaskInputs, TaskOutputs, Utils
from ..base_env.omix_env_task import BaseOmixEnvTask
from ..file.fasta_file import FastaFile
from ..file.gtf_file import GTFFile
from ..file.bam_file import BAMFile

@task_decorator("SalmonQuantStar")
class SalmonQuantStar(BaseOmixEnvTask):
    """
    SalmonQuantStar class. Represents a process that wraps Salmon quant program. This allowed to quantify, from a STAR mapping output file (.bam format), gene expressions.
    
    Configuration options
        * `threads`: Multi threading options: number of threads to use (min=6). [Default =  12].
        * `sequencing_strategy: To specify sequencing strategy (I = inward, O = outward, M = matching).
        * `library_protocol : To specify whether the protocol for the read library takes into account strand information [S = stranded U = unstranded]. [Default = U]"}
        See --> https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype.
    """

    input_specs = {
        'unsorted_transcriptome_star_bam_file': (BAMFile,),
        'fasta_genome_sequence': (FastaFile,),
        'gtf_annotation': (GTFFile,)            
    }
    output_specs = {
        'Salmon_quant_file': (File,)
    }
    config_specs = {
        "threads": IntParam(default_value=12, min_value=1, short_description="Number of threads [Default =  12] "),
        "sequencing_strategy": StrParam(allowed_values=["I", "O" , "M"], short_description="To specify sequencing strategy (I = inward, O = outward, M = matching). See https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype."),
        "library_protocol" : StrParam(allowed_values=["S", "U"], short_description="To specify whether the protocol for the read library takes into account strand information [S = stranded U = unstranded]. See https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype."),
    }
   
    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        result_file = File()
        result_file.path = self._get_output_file_path(params)
        return {"Salmon_quant_file": result_file}
    
    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        thread = params["threads"]
        genome_fasta = params["fasta_genome_sequence"]
        annotation = params["gtf_annotation"]  
        star_bam = params["unsorted_transcriptome_star_bam_file"]
        lib_strategy = params["library_protocol"]
        mapping_opt = params["sequencing_strategy"]
        experiment_opt = mapping_opt + "" + lib_strategy                        


        cmd = [
            " gffread -w transcripto.tmp.fa -g ",genome_fasta,
            annotation," ; cat transcripto.tmp.fa  | cut -d " " -f 1 > transcripto.tmp.2.fa ; rm transcripto.tmp.fa ;",
            "salmon quant  -p ",thread,
            " -t transcripto.tmp.2.fa -l ",experiment_opt,
            " -a ", star_bam,
            " -o ", self._get_output_file_path(params),
            " ; rm transcripto.tmp.2.fa ; "
        ] 

        return cmd

    def _get_output_file_path(self, params):
        return params["unsorted_transcriptome_star_bam_file"] + ".Salmon_counting"