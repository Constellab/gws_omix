# LICENSE
# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import File, Settings, BaseTestCase, TaskTester, GTest
from gws_omix import FastaFile
from gws_omix import GmapIndex
#from ..file.ec_list_file import ECListFile


class TestGmapIndex(BaseTestCase):
    
    async def test_GmapIndex(self):
        settings = Settings.retrieve()
        data_dir = settings.get_variable("gws_omix:testdata_dir")
        file = FastaFile()
        file.path = os.path.join(data_dir, "e_coli_K12.genome.fna.fasta")



        # run Gmap_build
        tester = TaskTester(
            inputs = {'uncompressed_genome_fasta_file': file},
            task_type = GmapIndex
        )
        outputs = await tester.run()
        f = outputs['gmap_index_file']
        result_content = os.path.join(f.path , "tmp.gmap.index.log")
        print(result_content.path)
        # # Get the expected file output       
        # expected_file_path = os.path.join(data_dir, "tmp.gmap.index.log")


        # with open(expected_file_path, 'r') as fp:
        #     expected_result_content = fp.read()
        #     print("----")
        #     print(expected_result_content)
        #     print("----")
        #     print(result_content)
        #     print("----")
        #     # Comparing results
        #     self.assertEqual( result_content , expected_result_content )
