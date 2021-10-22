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
        'fastq_dir': StrParam(default_value="", short_description="Location of the Fastq File") 
    }  
    input_specs = {

    }
    output_specs = {
        'fastq_file': (FastqFile,)
    }


    async def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        fastq = FastqFile()
        fastq.path = params['fastq_dir']
        return {"fastq_file": fastq}