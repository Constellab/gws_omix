# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import BaseTestCaseLight, File, Settings, TaskRunner
from gws_omix import BlastEC


# gws_omix/test_blast_ec
class TestBlastEC(BaseTestCaseLight):

    def test_BlastEC(self):
        data_dir = Settings.get_instance().get_variable("gws_omix:testdata_dir")
        file = File(os.path.join(data_dir, "fasta_subject.faa"))

        # Running BlastEC.py
        tester = TaskRunner(
            params={'taxonomy': 'fungi'},  # "uniprot_db_dir": os.path.join(testdata_dir, "uniprot_kb")
            inputs={'fasta_file': file},
            task_type=BlastEC
        )
        outputs = tester.run()
        blast_ec_file = outputs['filtered_blast_ec_file']
        result_content = blast_ec_file.read()

        # Get the expected file output
        expected_file_path = os.path.join(data_dir, "results.test.blast.csv")

        with open(expected_file_path, 'r') as fp:
            expected_result_content = fp.read()
            # Comparing results
            self.assertEqual(result_content, expected_result_content)
