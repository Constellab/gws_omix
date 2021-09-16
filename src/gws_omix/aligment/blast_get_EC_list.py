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


@task_decorator("BlastEcGetEcList")
class BlastEcGetEcList(BaseOmixEnvTask):
    """
    BlastEcGetEcList class. Get EC list File, from Blast_ec output File, for Digital Twins reconstructions.

    # blastec output file format (18 Columns) : qaccver saccver pident qcovs qcovhsp length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore alignment_type ec_number (NO header in file)
sp|P87228|SERA_SCHPO	sp|P87228|SERA_SCHPO	100.000	100	100	466	0	0	1	466	466	1	466	466	0.0	954	BEST_HIT	1.1.1.95; 1.1.1.399
sp|Q75DK1|SEC14_ASHGO	sp|Q75DK1|SEC14_ASHGO	100.000	100	100	308	0	0	1	308	308	1	308	308	0.0	637	BEST_HIT	NA
sp|Q75DK1|SEC14_ASHGO	sp|P24859|SEC14_KLULA	82.724	98	98	301	52	0	1	301	308	1	301	301	0.0	528	SECONDARY_HITS	NA
sp|Q75DK1|SEC14_ASHGO	sp|P53989|SEC14_CANGA	80.066	98	98	301	60	0	1	301	308	1	301	302	0.0	502	SECONDARY_HITS	NA
sp|A3LZ57|SEC16_PICST	sp|A3LZ57|SEC16_PICST	100.000	100	100	2212	0	0	1	2212	2212	1	2212	2212	0.0	4490	BEST_HIT	NA        

    """
    input_specs = {
        'EC_file': (BlastEcFile,)
    }
    output_specs = {
        'EC_list_output_file': (EcListFile,)
    }

   
    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        ec_fi  = params["EC_file"]
        result_file = File()
        result_file.path = self._create_ec_list_from_deepec_output_file(ec_fi)
        return {"EC_list_output_file": result_file} 

    def _create_ec_list_from_blast_output_file(self, blast_output_file):
        uniq_ec={}
        blast_ec = blast_output_file

        with open(blast_ec, 'r') as lines: 
            li=lines.readlines()
            for line in li:
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
