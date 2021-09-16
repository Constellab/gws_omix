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
from ..file.blast_ec_file import BlastEcFile
from ..file.ec_list_file import EcListFile

@task_decorator("BlastEcAndDeepEcGetEcList")
class BlastEcAndDeepEcGetEcList(BaseOmixEnvTask):

    input_specs = {
        'DeepEC_EC_file': (DeepEcFile,),
        'BLAST_EC_EC_file': (BlastEcFile,),        
    }
    output_specs = {
        'ec_list_merged': (EcListFile,)
    }
   
    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        deep_ec_fi  = params["DeepEC_EC_file"]
        blast_ec_fi  = params["BLAST_EC_EC_file"]
        result_file = File()
        result_file.path = self._create_ec_list_from_deepec_and_blast_ec_output_file(deep_ec_fi, blast_ec_fi)
        return {"EC_list_output_file": result_file} 


    def _create_ec_list_from_deepec_and_blast_ec_output_file(self, deepec_output_file, blast_output_file):
        uniq_ec={}
        deepec_ec = deepec_output_file
        blast_ec = blast_output_file

        with open(deepec_ec, 'r') as lines: 
            li_deepec=lines.readlines()
            for line in li_deepec:
                if re.match("^Query ID",line):
                    pass
                else:
                    li_split=line.split("\t")
                    ec_split=li_split[1].split(":")
                    uniq_ec[ec_split[1]]=1

        with open(blast_ec, 'r') as lines: 
            li_blast=lines.readlines()
            for line in li_blast:
                if re.match("^#",line):
                    pass
                else:
                    li_split=line.split("\t")
                    
                    if re.match("SECONDARY_HITS",str(li_split[16])):
                        pass                  
                    else:
                        ec_split=li_split[17].split("; ")
                        for ec in ec_split:
                            uniq_ec[ec]=1
        return uniq_ec
