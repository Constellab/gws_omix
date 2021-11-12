# LICENSE
# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import json
from gws_core import File, Settings, BaseTestCase, TaskTester, GTest, ViewTester, Folder, FileHelper
from gws_omix import SalmonIndex

class TestSalmonIndex(BaseTestCase):
    
    async def test_SalmonIndex(self):
#        SalmonIndex.uninstall()
        settings = Settings.retrieve()
        data_dir = settings.get_variable("gws_omix:testdata_dir")
        file = File()
        file.path = os.path.join(data_dir, "e_coli_K12.genome.fna.fasta")
        file_2 = File()
        file_2.path = os.path.join(data_dir, "test_ecoli.gtf")
        testdata_dir = settings.get_variable("gws_omix:testdata_dir")
  
        # Running 
        tester = TaskTester(
            params = {'threads': 2},
            inputs = {'uncompressed_genome_file': file, 'gtf_annotation': file_2 },
            task_type = SalmonIndex
        )
        outputs = await tester.run()

        f = outputs['salmon_index_folder']
        self.assertTrue(FileHelper.exists_on_os(f.path))
        file = File(path=os.path.join(f.path , "refInfo.json"))
        result_content = file.read()

        # Get the expected file output       
        expected_file_path = File(path=os.path.join(data_dir, "refInfo.json"))
        expected_result_content = expected_file_path.read()

        print("----")
        print(result_content)
        print("----")
        print(expected_result_content)
        self.assertEqual( result_content, expected_result_content  )
