# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, IntParam, StrParam, FloatParam, ConfigParams, TaskInputs, TaskOutputs, Utils, Task
from ..file.fasta_file import FastaFile


@task_decorator("FastaLoader")
class FastaLoader(Task):
    """
    FastaLoader Class. To load fasta file 

    """
    config_specs = {
        'fasta_dir': StrParam(default_value="", short_description="Location of the Fasta File") 
    }  
    input_specs = {

    }
    output_specs = {
        'fasta_file': (FastaFile,)
    }


    async def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        fasta = FastaFile()
        fasta.path = params['fasta_dir']
        return {"fasta_file": fasta}