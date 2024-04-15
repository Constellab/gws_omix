

import json
import os

from gws_core import BaseTestCase, Settings, TaskRunner, ViewTester
from gws_omix import BlastECFile, BlastECListExtractor


# test_blast_ec_list_extractor
class TestBlastEC(BaseTestCase):

    def test_BlastECListExtractor(self):
        settings = Settings.retrieve()
        data_dir = settings.get_variable("gws_omix:testdata_dir")
        file = BlastECFile()
        file.path = os.path.join(data_dir, "results.test.blast.csv")

        # run BlastECListExtractor
        tester = TaskRunner(
            inputs={'blastec_file': file},
            task_type=BlastECListExtractor
        )
        outputs = tester.run()
        f = outputs['ec_list_file']
        result_content = f.read()

        # Get the expected file output
        expected_file_path = os.path.join(
            data_dir, "results.test.blast.ec_list.txt")

        with open(expected_file_path, 'r') as fp:
            expected_result_content = fp.read()
            print("----")
            print(result_content)
            print("----")
            print(expected_result_content)
            print("----")
            # Comparing results
            self.assertEqual(result_content, expected_result_content)

        text_view = f.view_as_raw_text()
        tester = ViewTester(view=text_view)
        text_view_dict = tester.to_dict(params={"page": 1, "page_size": 1000})
        print(json.dumps(text_view_dict, indent=2))
        self.assertEqual(text_view_dict["total_number_of_pages"], 1)
