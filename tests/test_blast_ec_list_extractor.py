# LICENSE
# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import json

from gws_core import File, Settings, BaseTestCase, TaskRunner, GTest, ViewTester
from gws_omix import BlastECListExtractor
from gws_omix import BlastECFile
#from ..file.ec_list_file import ECListFile


class TestBlastEC(BaseTestCase):
    
    async def test_BlastECListExtractor(self):
        settings = Settings.retrieve()
        data_dir = settings.get_variable("gws_omix:testdata_dir")
        file = BlastECFile()
        file.path = os.path.join(data_dir, "results.test.blast.csv")

        # run BlastECListExtractor
        tester = TaskRunner(
            inputs = {'blastec_file': file},
            task_type = BlastECListExtractor
        )
        outputs = await tester.run()
        f = outputs['ec_list_file']
        result_content = f.read()

        # Get the expected file output       
        expected_file_path = os.path.join(data_dir, "results.test.blast.ec_list.txt")

        with open(expected_file_path, 'r') as fp:
            expected_result_content = fp.read()
            print("----")
            print(result_content)
            print("----")
            print(expected_result_content)
            print("----")
            # Comparing results
            self.assertEqual( result_content, expected_result_content  )

        text_view = f.view_as_raw_text()
        tester = ViewTester(view=text_view)
        text_view_dict = tester.to_dict(params={"page":1, "page_size":1000})
        print(json.dumps(text_view_dict, indent=2))
        self.assertEqual( text_view_dict["total_number_of_pages"], 1 )