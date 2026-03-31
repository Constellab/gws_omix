import os

from gws_core import BaseTestCase, File, Folder, Settings, TaskRunner

from gws_omix import FastqFolder
from gws_omix.rna_seq.mapping_transcriptome.salmon import SalmonQuant
from gws_omix.rna_seq.mapping_transcriptome.salmon_index import SalmonIndex

TESTDATA = os.path.join(Settings.get_instance().get_variable("gws_omix", "testdata_dir"), "rna_seq")


class TestSalmon(BaseTestCase):

    def test_salmon_index(self):
        """SalmonIndex builds an index from mini transcript FASTA."""
        self.print("Test SalmonIndex: mini_transcripts.fa")

        output: Folder = TaskRunner(
            task_type=SalmonIndex,
            inputs={"transcript_fasta": File(os.path.join(TESTDATA, "mini_transcripts.fa"))},
            params={"threads": 2},
        ).run()["output"]

        self.assertIsNotNone(output)
        index_files = os.listdir(output.path)
        self.print(f"Salmon index files: {index_files}")
        self.assertGreater(len(index_files), 0, msg="Expected Salmon index files in output")

    def test_salmon_quant_se(self):
        """Chain SalmonIndex → SalmonQuant SE: produces quant.sf per sample."""
        self.print("Test SalmonQuant SE: build index then quantify mini_SE.fastq.gz")

        idx: Folder = TaskRunner(
            task_type=SalmonIndex,
            inputs={"transcript_fasta": File(os.path.join(TESTDATA, "mini_transcripts.fa"))},
            params={"threads": 2},
        ).run()["output"]

        quant_folder: Folder = TaskRunner(
            task_type=SalmonQuant,
            inputs={
                "fastq_folder": FastqFolder(TESTDATA),
                "metadata": File(os.path.join(TESTDATA, "metadata_se.tsv")),
                "salmon_index": idx,
            },
            params={"threads": 2},
        ).run()["output"]

        self.assertIsNotNone(quant_folder)
        quant_sf_files = []
        for entry in os.listdir(quant_folder.path):
            qsf = os.path.join(quant_folder.path, entry, "quant.sf")
            if os.path.isfile(qsf):
                quant_sf_files.append(qsf)
        self.print(f"quant.sf files: {quant_sf_files}")
        self.assertGreater(len(quant_sf_files), 0, msg="Expected quant.sf in Salmon output")

    def test_salmon_quant_pe(self):
        """Chain SalmonIndex → SalmonQuant PE: produces quant.sf for paired-end sample."""
        self.print("Test SalmonQuant PE: build index then quantify mini_R1/R2.fastq.gz")

        idx: Folder = TaskRunner(
            task_type=SalmonIndex,
            inputs={"transcript_fasta": File(os.path.join(TESTDATA, "mini_transcripts.fa"))},
            params={"threads": 2},
        ).run()["output"]

        quant_folder: Folder = TaskRunner(
            task_type=SalmonQuant,
            inputs={
                "fastq_folder": FastqFolder(TESTDATA),
                "metadata": File(os.path.join(TESTDATA, "metadata_pe.tsv")),
                "salmon_index": idx,
            },
            params={"threads": 2},
        ).run()["output"]

        self.assertIsNotNone(quant_folder)
        quant_sf_files = []
        for entry in os.listdir(quant_folder.path):
            qsf = os.path.join(quant_folder.path, entry, "quant.sf")
            if os.path.isfile(qsf):
                quant_sf_files.append(qsf)
        self.print(f"PE quant.sf files: {quant_sf_files}")
        self.assertGreater(len(quant_sf_files), 0, msg="Expected quant.sf in Salmon PE output")
