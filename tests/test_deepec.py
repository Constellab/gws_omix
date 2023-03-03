# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os

from gws_core import (BaseTestCase, ConfigParams, Settings, TaskRunner,
                      ViewTester)
from gws_omix import DeepEC, FastaFile

# from ..file.ec_list_file import ECListFile


class TestDeepEc(BaseTestCase):

    def test_DeepEc(self):
        settings = Settings.retrieve()
        data_dir = settings.get_variable("gws_omix:testdata_dir")
        file = FastaFile()
        file.path = os.path.join(data_dir, "fasta.test.deepec.fasta")
        print(file.path)
        # run DeepEC
        tester = TaskRunner(
            inputs={'fasta_file': file},
            task_type=DeepEC
        )
        outputs = tester.run()
        f = outputs['deepec_file']
        result_content = f.read()

        # Get the expected file output
        expected_file_path = os.path.join(data_dir, "log_files/DeepEC_Result_DL.txt")

        with open(expected_file_path, "r", encoding="utf-8") as fp:
            expected_result_content = fp.read()
            print("----")
            print(result_content)
            print("----")
            print(expected_result_content)
            print("----")
            # Comparing results
            self.assertEqual(result_content, expected_result_content)

        text_view = f.view_as_raw_text(ConfigParams())
        tester = ViewTester(view=text_view)
        text_view_dict = tester.to_dict(params={"page": 1, "page_size": 1000})
        print(json.dumps(text_view_dict, indent=2))
        self.assertEqual(text_view_dict["data"]["total_number_of_pages"], 1)
        self.assertEqual(text_view_dict["data"]["total_number_of_items"], 304)
