import os

import pandas as pd
from gws_core import BaseTestCase, File, Folder, Settings, TaskRunner

from gws_omix.rna_seq.mapping_transcriptome.salmon_merge_matrix import SalmonMergeMatrix

TESTDATA = os.path.join(Settings.get_instance().get_variable("gws_omix", "testdata_dir"), "rna_seq")


class TestSalmonMergeMatrix(BaseTestCase):

    def test_merge_synthetic_quant_sf(self):
        """SalmonMergeMatrix merges 2 synthetic quant.sf into a transcript count matrix."""
        self.print("Test SalmonMergeMatrix: synthetic salmon_quant/ (sample1 + sample2)")

        matrix_file: File = TaskRunner(
            task_type=SalmonMergeMatrix,
            inputs={"salmon_quant_folder": Folder(os.path.join(TESTDATA, "salmon_quant"))},
            params={},
        ).run()["matrix"]

        self.assertIsNotNone(matrix_file)
        self.assertTrue(os.path.exists(matrix_file.path))

        df = pd.read_csv(matrix_file.path)
        self.print(f"Merged matrix shape: {df.shape}")
        self.print(df.to_string())

        self.assertIn("transcript_id", df.columns, msg="Expected 'transcript_id' column")
        # 2 sample columns + transcript_id
        self.assertGreaterEqual(df.shape[1], 3, msg="Expected at least 3 columns")
        # 3 transcripts in synthetic quant.sf
        self.assertEqual(df.shape[0], 3, msg="Expected 3 transcripts in merged matrix")
