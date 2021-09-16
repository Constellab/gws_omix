# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, StrParam, ConfigParams, TaskInputs, TaskOutputs, Utils
from ..base.omix_env_task import BaseOmixEnvTask
from ..file.deepec_file import DeepEcFile
from ..file.ec_list_file import EcListFile


@task_decorator("DeepEcGetEcList")
class DeepEcGetEcList(BaseOmixEnvTask):
    """
    DeepEcGetEcList class. Get EC list File, from DeepEC output File, for Digital Twins reconstructions.

    # deepec output file format (2 Columns): Query ID and Predicted ec number (header in file)
    Query ID   Predicted EC number
    Mycgr3P87220	EC:2.7.1.24
    Mycgr3P42503	EC:3.2.1.3
    Mycgr3P66454	EC:3.1.1.3
    Mycgr3P65205	EC:1.2.1.3
        
    """

    input_specs = {
        'EC_file': (DeepEcFile,)
    }
    output_specs = {
        'EC_list_output_file': (EcListFile,)
    }

   
    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        result_file = File()
        result_file.path = self._create_ec_list_from_deepec_output_file(EC_file)
        return {"EC_list_output_file": result_file} 

# deepec
#Query ID   Predicted EC number
#Mycgr3P87220	EC:2.7.1.24
#Mycgr3P42503	EC:3.2.1.3
#Mycgr3P66454	EC:3.1.1.3
#Mycgr3P65205	EC:1.2.1.3
    def _create_ec_list_from_deepec_output_file(self, deepec_output_file):
        uniq_ec={}
        deepec_ec = deepec_output_file

        with open(deepec_ec, 'r') as lines: 
            li=lines.readlines()
            for line in li:
                if re.match("^Query ID",line):
                    pass
                else:
                    li_split=line.split("\t")
                    ec_split=li_split[1].split(":")
                    uniq_ec[ec_split[1]]=1
        return uniq_ec
