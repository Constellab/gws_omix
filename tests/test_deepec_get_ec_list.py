# LICENSE
# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import json

from gws_core import File, Settings, BaseTestCase, TaskTester, GTest
from gws_omix import DeepEcGetEcList
from gws_omix import DeepECFile
#from ..file.ec_list_file import ECListFile


class TestDeepEC(BaseTestCase):
    
    async def test_BlastEcGetEcList(self):
        settings = Settings.retrieve()
        data_dir = settings.get_variable("gws_omix:testdata_dir")
        file = DeepECFile()
        file.path = os.path.join(data_dir, "DeepEC_Result.txt")



        # run BlastEcGetEcList
        tester = TaskTester(
            inputs = {'deepec_file': file},
            task_type = DeepEcGetEcList
        )
        outputs = await tester.run()
        f = outputs['ec_list_file']
        result_content = f.read()

        # Get the expected file output       
        expected_file_path = os.path.join(data_dir, "results.test.deepec.ec_list.txt")

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
        text_view_dict = text_view.to_dict(page=1, page_size=500)
        print(json.dumps(text_view_dict, indent=2))
        self.assertEqual( text_view_dict["total_number_of_pages"], 1 )
