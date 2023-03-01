# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com
import os
import re

from gws_core import (InputSpec, OutputSpec, TaskInputs,
                      TaskOutputs, task_decorator, ConfigParams, InputSpec, OutputSpec, InputSpecs, OutputSpecs)

from ..base_env.omix_env_task import BaseOmixEnvTask
from ..file.blast_ec_file import BlastECFile
from ..file.ec_list_file import ECListFile


@task_decorator("BlastECListExtractor",
                short_description="Extracts an EC list from a BlastECFile")
class BlastECListExtractor(BaseOmixEnvTask):
    """
    BlastECListExtractor class.

    Extracts an EC list from a `BlastECFile` for Digital Twins reconstructions.
    Keeps only the best hit for each previously analysed sequence.

    # BlastECListExtractor output file format (18 Columns) : qaccver saccver pident qcovs qcovhsp length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore alignment_type ec_number (NO header in file)
    sp|P87228|SERA_SCHPO	sp|P87228|SERA_SCHPO	100.000	100	100	466	0	0	1	466	466	1	466	466	0.0	954	BEST_HIT	1.1.1.95; 1.1.1.399
    sp|Q75DK1|SEC14_ASHGO	sp|Q75DK1|SEC14_ASHGO	100.000	100	100	308	0	0	1	308	308	1	308	308	0.0	637	BEST_HIT	NA
    sp|Q75DK1|SEC14_ASHGO	sp|P24859|SEC14_KLULA	82.724	98	98	301	52	0	1	301	308	1	301	301	0.0	528	SECONDARY_HITS	NA
    sp|Q75DK1|SEC14_ASHGO	sp|P53989|SEC14_CANGA	80.066	98	98	301	60	0	1	301	308	1	301	302	0.0	502	SECONDARY_HITS	NA
    sp|A3LZ57|SEC16_PICST	sp|A3LZ57|SEC16_PICST	100.000	100	100	2212	0	0	1	2212	2212	1	2212	2212	0.0	4490	BEST_HIT	NA

    """
    input_specs: InputSpecs = {
        'blastec_file': InputSpec(BlastECFile, human_name="", short_description="")
    }
    output_specs: OutputSpecs = {
        'ec_list_file': OutputSpec(ECListFile, human_name="", short_description="")
    }

    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        ec_fi = inputs["blastec_file"]
        result_file = ECListFile()
        result_file.path = self._create_ec_list_from_blast_output_file(ec_fi)
        return {"ec_list_file": result_file}

    def _create_ec_list_from_blast_output_file(self, blast_output_file):
        uniq_ec = {}
        blast_output_file_name = os.path.basename(blast_output_file.path)
        ec_list_path = self._get_output_file_path(blast_output_file_name)
        with open(ec_list_path, 'w+') as ec_list:
            with open(blast_output_file.path, 'r') as lines:
                li = lines.readlines()
                for line in li:
                    if re.match("^#", line):
                        pass
                    else:
                        li_split = line.split("\t")
                        if re.match("SECONDARY_HITS", str(li_split[16])):
                            pass
                        else:
                            if re.search(';', str(li_split[17])):
                                ec_split = re.split(';\s', li_split[17])
                                #ec_list.write("\n".join(str(ec) for ec in ec_split))
                                ec_list.write("\n".join(li_split[0]+"\t"+str(ec) for ec in ec_split))
                            elif re.match("NA", str(li_split[17])):
                                pass
                            else:
                                #ec_split = li_split[17]
                                ec_split = li_split[0]+"\t"+li_split[17]
                                ec_list.write(ec_split)
#                                uniq_ec[ec]=1
                # for ec in uniq_ec:
                #     ec_list.write( ec )
        return ec_list_path

    def _get_output_file_path(self, blast_output_file_name):
        return os.path.join(
            self.working_dir,
            blast_output_file_name + ".blast.ec_list.txt"
        )
