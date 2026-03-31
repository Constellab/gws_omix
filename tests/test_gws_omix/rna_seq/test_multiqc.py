import os

from gws_core import BaseTestCase, File, Folder, Settings, TaskRunner

from gws_omix import FastqFolder
from gws_omix.rna_seq.multiqc.multiqc import MultiQc
from gws_omix.rna_seq.quality_check.fastq_init import FastqcInit

TESTDATA = os.path.join(Settings.get_instance().get_variable("gws_omix", "testdata_dir"), "rna_seq")


class TestMultiQc(BaseTestCase):

    def test_multiqc_from_fastqc_reports(self):
        """Chain FastQC → MultiQC: MultiQC aggregates FastQC reports into a single HTML."""
        self.print("Test MultiQC: FastQC SE → MultiQC")

        # Step 1: run FastQC to produce reports
        fastqc_output: Folder = TaskRunner(
            task_type=FastqcInit,
            inputs={
                "fastq_folder": FastqFolder(TESTDATA),
                "metadata": File(os.path.join(TESTDATA, "metadata_se.tsv")),
            },
            params={"threads": 2},
        ).run()["output"]
        self.print(f"FastQC reports dir: {os.listdir(fastqc_output.path)}")

        # Step 2: run MultiQC on FastQC output
        output: Folder = TaskRunner(
            task_type=MultiQc,
            inputs={"fastqc_reports_folder": fastqc_output},
            params={},
        ).run()["output"]

        self.assertIsNotNone(output)
        files = os.listdir(output.path)
        self.print(f"MultiQC output files: {files}")
        html_files = [f for f in files if f.endswith(".html")]
        self.assertGreater(len(html_files), 0, msg="Expected MultiQC combined HTML report")
