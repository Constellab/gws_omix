import os

from gws_core import BaseTestCase, File, Folder, Settings, TaskRunner

from gws_omix import FastqFolder
from gws_omix.rna_seq.quality_check.fastq_init import FastqcInit

TESTDATA = os.path.join(Settings.get_instance().get_variable("gws_omix", "testdata_dir"), "rna_seq")


class TestFastqcInit(BaseTestCase):

    def test_fastqc_single_end(self):
        """FastQC runs on SE FASTQ and produces HTML reports."""
        self.print("Test FastQC SE: mini_SE.fastq.gz")
        runner = TaskRunner(
            task_type=FastqcInit,
            inputs={
                "fastq_folder": FastqFolder(TESTDATA),
                "metadata": File(os.path.join(TESTDATA, "metadata_se.tsv")),
            },
            params={"threads": 2},
        )
        output: Folder = runner.run()["output"]
        self.assertIsNotNone(output)
        files = os.listdir(output.path)
        self.print(f"FastQC output files: {files}")
        html_files = [f for f in files if f.endswith(".html")]
        self.assertGreater(len(html_files), 0, msg="Expected at least one FastQC HTML report")

    def test_fastqc_paired_end(self):
        """FastQC runs on PE FASTQ (R1+R2) and produces 2 reports."""
        self.print("Test FastQC PE: mini_R1 + mini_R2")
        runner = TaskRunner(
            task_type=FastqcInit,
            inputs={
                "fastq_folder": FastqFolder(TESTDATA),
                "metadata": File(os.path.join(TESTDATA, "metadata_pe.tsv")),
            },
            params={"threads": 2},
        )
        output: Folder = runner.run()["output"]
        self.assertIsNotNone(output)
        files = os.listdir(output.path)
        html_files = [f for f in files if f.endswith(".html")]
        self.print(f"FastQC PE HTML reports: {html_files}")
        self.assertGreaterEqual(len(html_files), 2, msg="Expected reports for both R1 and R2")
