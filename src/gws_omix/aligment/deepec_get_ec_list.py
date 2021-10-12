# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, StrParam, ConfigParams, TaskInputs, TaskOutputs, Utils
from ..base.omix_env_task import BaseOmixEnvTask
from ..file.deepec_file import DeepECFile
from ..file.ec_list_file import ECListFile

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
        'deepec_file': (DeepECFile,)
    }
    output_specs = {
        'ec_list_file': (ECListFile,)
    }
   
    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        ec_fi  = inputs["deepec_file"]
        result_file = ECListFile()
        result_file.path = self._create_ec_list_from_deepec_output_file(ec_fi)
        return {"ec_list_file": result_file}

    def _create_ec_list_from_deepec_output_file(self, deepec_output_file):
        uniq_ec={}
        deepec_output_file_name  = os.path.basename(deepec_output_file.path)
        ec_list_path = self._get_output_file_path(deepec_output_file_name)
        with open(ec_list_path, 'w+') as ec_list: 
            with open(deepec_output_file.path, 'r') as lines: 
                li=lines.readlines()
                for line in li:
                    if re.match("^Query ID",line):
                        pass
                    else:
                        li_split=line.split("\t")
                        ec_split=li_split[1].split(":")
                        uniq_ec[ec_split[1]]=1
                for ec in uniq_ec:
                    ec_list.write( ec )
        return ec_list_path

    def _get_output_file_path(self, deepec_output_file_name) :
        return os.path.join(
            self.working_dir, 
            deepec_output_file_name + ".deepec.ec_list.txt"
        )



#     def _create_ec_list_from_blast_output_file(self, blast_output_file):
#         uniq_ec={}
#         blast_output_file_name  = os.path.basename(blast_output_file.path)
#         ec_list_path = self._get_output_file_path(blast_output_file_name)
#         with open(ec_list_path, 'w+') as ec_list: 
#             with open(blast_output_file.path , 'r') as lines: 
#                 li=lines.readlines()
#                 for line in li:
#                     if re.match("^#",line):
#                         pass
#                     else:
#                         li_split=line.split("\t")                    
#                         if re.match("SECONDARY_HITS",str(li_split[16])):
#                             pass                  
#                         else:
#                             if re.search(';',str(li_split[17])):
#                                 ec_split=re.split(';\s',li_split[17])
#                                 ec_list.write("\n".join(str(ec) for ec in ec_split))
#                             elif re.match("NA",str(li_split[17])):
#                                 pass
#                             else:
#                                 ec_split=li_split[17]
#                                 ec_list.write( ec_split)
# #                                uniq_ec[ec]=1
#                 # for ec in uniq_ec:
#                 #     ec_list.write( ec )
#         return ec_list_path

#     def _get_output_file_path(self, blast_output_file_name) :
#         return os.path.join(
#             self.working_dir, 
#             blast_output_file_name + ".blast.ec_list.txt"
#         )
