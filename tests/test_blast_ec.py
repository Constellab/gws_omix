# LICENSE
# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import json
from gws_core import File, Settings, BaseTestCase, TaskTester, GTest, ViewTester
from gws_omix import BlastEC

class TestBlastEC(BaseTestCase):
    
    async def test_BlastEC(self):
        BlastEC.uninstall()
        settings = Settings.retrieve()
        data_dir = settings.get_variable("gws_omix:testdata_dir")
        file = File()
        file.path = os.path.join(data_dir, "fasta_subject.faa")

        settings = Settings.retrieve()
        testdata_dir = settings.get_variable("gws_omix:testdata_dir")

        # Running BlastEC.py
        tester = TaskTester(
            params = {'taxonomy': 'fungi', "uniprot_db_dir": os.path.join(testdata_dir, "uniprot_kb")},
            inputs = {'fasta_file': file},
            task_type = BlastEC
        )
        outputs = await tester.run()
        blast_ec_file = outputs['filtered_blast_ec_file']
        result_content = blast_ec_file.read()


        # Get the expected file output       
        expected_file_path = os.path.join(data_dir, "results.test.blast.csv")

        with open(expected_file_path, 'r') as fp:
            expected_result_content = fp.read()
            print("----")
            print(result_content)
            print("----")
            print(expected_result_content)
            print("----")
            # Comparing results
            self.assertEqual( result_content, expected_result_content  )


        table_view = blast_ec_file.view_head_as_table(params=None)
        tester = ViewTester(view=table_view)
        dic_ = tester.to_dict(params={})
        print(json.dumps(dic_, indent=2))
        self.assertEqual( dic_["total_number_of_rows"], 24 )
