# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import re

from gws_core import BaseTestCase, File, FileHelper, Settings, TaskRunner
from gws_omix import SalmonIndex


class TestSalmonIndex(BaseTestCase):

    def test_SalmonIndex(self):
        # SalmonIndex.uninstall()
        settings = Settings.retrieve()
        data_dir = settings.get_variable("gws_omix:testdata_dir")
        file = File()
        file.path = os.path.join(data_dir, "e_coli_K12.genome.fna.fasta")
        file_2 = File()
        file_2.path = os.path.join(data_dir, "test_ecoli.gtf")
        testdata_dir = settings.get_variable("gws_omix:testdata_dir")

        # Running
        tester = TaskRunner(
            params={'threads': 2},
            inputs={'uncompressed_genome_file': file, 'gtf_annotation': file_2},
            task_type=SalmonIndex
        )
        outputs = tester.run()

        f = outputs['salmon_index_result']
        self.assertTrue(FileHelper.exists_on_os(f.path))
        file = File(path=os.path.join(f.path, "pre_indexing.log"))
        result_content = file.read()
        pattern = re.compile("done building index")
        result_file_ok = ""
        for line in result_content:
            for match in re.finditer(pattern, line):
                result_file_ok = "OK"

        # Get the expected file output
        expected_file_path = File(path=os.path.join(
            data_dir, "e_coli_K12.genome.fna.fasta.salmon_index/pre_indexing.log"))
        expected_result_content = expected_file_path.read()
        expected_result_file_ok = ""
        for line_2 in expected_result_content:
            for match in re.finditer(pattern, line_2):
                expected_result_file_ok = "OK"

        print("----")
        print(result_content)
        print("----")
        print(expected_result_content)
        self.assertEqual(result_file_ok, expected_result_file_ok)
