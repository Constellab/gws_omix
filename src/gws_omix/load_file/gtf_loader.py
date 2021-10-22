# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, IntParam, StrParam, FloatParam, ConfigParams, TaskInputs, TaskOutputs, Utils, Task
from ..file.gtf_file import GTFFile


@task_decorator("GtfLoader")
class GtfLoader(Task):
    """
    GtfLoader Class. To load GTF annotation file 

    """
    config_specs = {
        'gtf_dir': StrParam(default_value="", short_description="Location of the gtf annotation File") 
    }  
    input_specs = {

    }
    output_specs = {
        'gtf_file': (GTFFile,)
    }


    async def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        gtf = GTFFile()
        gtf.path = params['gtf_dir']
        return {"gtf_file": gtf}