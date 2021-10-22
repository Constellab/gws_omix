# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, IntParam, StrParam, FloatParam, ConfigParams, TaskInputs, TaskOutputs, Utils, Task
from ..file.gff3_file import GFF3File


@task_decorator("Gff3Loader")
class Gff3Loader(Task):
    """
    Gff3Loader Class. To load GFF3 annotation file 

    """
    config_specs = {
        'gff3_dir': StrParam(default_value="", short_description="Location of the gff3 File") 
    }  
    input_specs = {

    }
    output_specs = {
        'gff3_file': (GFF3File,)
    }


    async def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        gff3 = GFF3File()
        gff3.path = params['gff3_dir']
        return {"gff3_file": gff3}