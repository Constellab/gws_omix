# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, IntParam, StrParam, ConfigParams, TaskInputs, TaskOutputs, Utils
from ..base.omix_env_task import BaseOmixEnvTask
from ..file.salmon_tpm_quantmerge_output_file import SalmonTpmQuantmergeOutputFile
from ..file.salmon_reads_quantmerge_output_file import SalmonReadsQuantmergeOutputFile



@task_decorator("SalmonQuantMerge")
class SalmonQuantMerge(BaseOmixEnvTask):
    """
    Salmon quantmerge class. Represents a process that wraps Salmon quantmerge program.

    Allowed to obtain tabular files which includes a merging of previous salmon quant output files (raw reads count or normalized TPM)

    # sample.raw_reads_output.tsv
    Name    sample_1    sample_2    sample_3
    gene_1  <Raw reads count>   <Raw reads count>   <Raw reads count>
    gene_2  <Raw reads count>   <Raw reads count>   <Raw reads count>
    gene_3  <Raw reads count>   <Raw reads count>   <Raw reads count>
    
    # sample.tpm_output.tsv
    Name    sample_1    sample_2    sample_3
    gene_1  <Raw reads count>   <Raw reads count>   <Raw reads count>
    gene_2  <Raw reads count>   <Raw reads count>   <Raw reads count>
    gene_3  <Raw reads count>   <Raw reads count>   <Raw reads count>   

    """

    output_specs = { 
        'Salmon_quant_TPM_file': (SalmonTpmQuantmergeOutputFile,),
        'Salmon_quant_RAW_reads_file': (SalmonReadsQuantmergeOutputFile,)
     }
    config_specs = {
        "experiment_name": StrParam(default_value="Current_experiment",  short_description=" Output file names "),
        "salmon_quant_output_directories" : StrParam(short_description=" The directory list (output directory from salmon quant) can be written in two ways : ./sample_1 ./sample_2 ... ./sample_n |OR| ./sample_* (for instance, if the first word(s) are shared for ervery samples) ")
    }
   
    def gather_output(self, stdout: str=None):
        self.output["Salmon_quant_TPM_file"] = File( path=self.output_file_TPM )
        self.output["Salmon_quant_RAW_reads_file"] = File( path=self.output_file_RAW_count )
    
    def build_command(self) -> list:
        brick_dir = self.get_brick_dir("omix")
        bin_file = os.path.join(brick_dir, "bin", "Salmon")   

        exp_name = self.get_param("experiment_name")
        self.output_file_TPM = os.path.join(exp_name,".tpm_output.tsv")
        self.output_file_RAW_count = os.path.join(exp_name,".raw_reads_output.tsv")       
        cmd = [
            bin_file,
            " quantmerge --quants *.RNAseq_mapping.genome.Salmon_counting --column tpm -o ",self.output_file_TPM,
            " ; ",
            bin_file,
            " quantmerge --quants *.RNAseq_mapping.genome.Salmon_counting --column numreads -o ",self.output_file_RAW_count,
            " ; "
        ]      

        return cmd