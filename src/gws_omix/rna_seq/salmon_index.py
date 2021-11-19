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
from ..file.salmon_index_result_folder import SalmonIndexResultFolder

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
        'salmon_index_result': (SalmonIndexResultFolder,)
    }
    config_specs = {
        "threads": IntParam(default_value=2, min_value=1, short_description="Number of threads [Default =  2] ")
    }
    
    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        path = self._get_output_file_path(inputs)

        print("xxxxx")
        print(path)
        print(os.listdir(os.path.join(path)))
        print("---")
        with open(os.path.join(path, "ref_indexing.log")) as fp:
            print(fp.read())
        print("---")

        result_file = SalmonIndexResultFolder(path=path)
        return {"salmon_index_result": result_file} 
    
    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        thread = params["threads"]
        annotation = inputs["gtf_annotation"]
        genome_fasta = inputs["uncompressed_genome_file"]
        script_file_dir = os.path.dirname(os.path.realpath(__file__))

        cmd = [
            "bash",
            os.path.join(script_file_dir, "./sh/salmon_index_cmd.sh"),
            genome_fasta.path,
            annotation.path,
            thread,
            self._get_fasta_file_name(inputs) 
        ]
        return cmd
    
    def _get_fasta_file_name(self, inputs):
        genome_fasta = inputs["uncompressed_genome_file"]
        return os.path.basename(genome_fasta.path)

    def _get_output_file_path(self, inputs):
        return os.path.join(
            self.working_dir, 
            self._get_fasta_file_name(inputs) + ".salmon_index"
        )
