# LICENSE
# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import json
from gws_core import File, Settings, BaseTestCase, TaskTester, GTest
from gws_omix import GmapAlignGFF3
from gws_omix import GFF3File
from gws_omix import FastaFile

class TestGmapGFF3(BaseTestCase):
    
    async def test_GmapGFF3(self):
        settings = Settings.retrieve()
        data_dir = settings.get_variable("gws_omix:testdata_dir")
        file_1 = FastaFile()
        file_1.path = os.path.join(data_dir, "e_coli_K12.genome.fna.fasta")
        file_2 = FastaFile()
        file_2.path = os.path.join(data_dir, "test.gmap.fasta.fna")

#"cross-species":"Yes", "alt-start-codons":"Yes", "fulllength":"Yes"
        # Running BlastEC.py
        tester = TaskTester(
            params = {"threads":2, "min-identity":0.0, "min-trimmed-coverage":0.0, "max-hit-number":1 },
            inputs = {'uncompressed_genome_fasta_file': file_1, 'cdna_or_cds_fasta_file': file_2},
            task_type = GmapAlignGFF3
        )

        outputs = await tester.run()
        gff3_file = outputs['gmap_gff3_file']
        result_content = gff3_file.read()


        # Get the expected file output       
        expected_file_path = os.path.join(data_dir, "test.gmap.gff3")

        with open(expected_file_path, 'r') as fp:
            expected_result_content = fp.read()
            print("----")
            print(result_content)
            print("----")
            print(expected_result_content)
            print("----")
            # Comparing results
            self.assertEqual( result_content, expected_result_content  )


        text_view = gff3_file.view_head_as_raw_text()
        text_view_dict = text_view.to_dict(page=1, page_size=500)
        print(json.dumps(text_view_dict, indent=2))
        self.assertEqual( text_view_dict["total_number_of_pages"], 3 )


