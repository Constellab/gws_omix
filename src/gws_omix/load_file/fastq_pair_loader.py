# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, IntParam, StrParam, FloatParam, ConfigParams, TaskInputs, TaskOutputs, Utils, Task
from ..file.fastq_file import FastqFile


@task_decorator("FastqLoader")
class FastqLoader(Task):
    """
    FastqLoader Class. To load fasta file 

    """
    config_specs = {
        'fastq_2_dir': StrParam(default_value="", short_description="Location of the Fastq File forward (_1)"),
        'fastq_1_dir': StrParam(default_value="", short_description="Location of the Fastq File reverse (_2)")
    }  
    input_specs = {

    }
    output_specs = {
        'fastq_1_file': (FastqFile,),
        'fastq_2_file': (FastqFile,)
    }


    async def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        fastq = FastqFile()
        fastq_1.path = params['fastq_1_dir']
        fastq_2.path = params['fastq_2_dir']
        return {"fastq_1_file": fastq_1, "fastq_2_file": fastq_2}