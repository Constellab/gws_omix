import os

from gws_core import BaseTestCase, File, Settings, TaskRunner

from gws_omix import FastqFolder
from gws_omix.rna_seq.trimming.fastp import Fastp

TESTDATA = os.path.join(Settings.get_instance().get_variable("gws_omix", "testdata_dir"), "rna_seq")


class TestFastp(BaseTestCase):

    def test_fastp_single_end(self):
        """Fastp trims SE FASTQ and returns FastqFolder with trimmed files."""
        self.print("Test Fastp SE: mini_SE.fastq.gz")
        runner = TaskRunner(
            task_type=Fastp,
            inputs={
                "fastq_folder": FastqFolder(TESTDATA),
                "metadata": File(os.path.join(TESTDATA, "metadata_se.tsv")),
            },
            params={
                "threads": 2,
                "Forward_separator": "",
                "Reverse_separator": "",
                "5_prime_hard_trimming_read_size": 0,
            },
        )
        output = runner.run()["output"]
        self.assertIsNotNone(output)
        files = os.listdir(output.path)
        self.print(f"Fastp SE output files: {files}")
        trimmed = [f for f in files if f.endswith(".fastq.gz")]
        self.assertGreater(len(trimmed), 0, msg="Expected trimmed FASTQ in output")

    def test_fastp_paired_end(self):
        """Fastp trims PE FASTQ (R1+R2) and produces trimmed FastqFolder."""
        self.print("Test Fastp PE: mini_R1.fastq.gz + mini_R2.fastq.gz")
        runner = TaskRunner(
            task_type=Fastp,
            inputs={
                "fastq_folder": FastqFolder(TESTDATA),
                "metadata": File(os.path.join(TESTDATA, "metadata_pe.tsv")),
            },
            params={
                "threads": 2,
                "Forward_separator": "R1",
                "Reverse_separator": "R2",
                "5_prime_hard_trimming_read_size": 0,
            },
        )
        output = runner.run()["output"]
        self.assertIsNotNone(output)
        files = os.listdir(output.path)
        self.print(f"Fastp PE output files: {files}")
        trimmed = [f for f in files if f.endswith(".fastq.gz")]
        self.assertGreater(len(trimmed), 0, msg="Expected trimmed FASTQ for PE")
