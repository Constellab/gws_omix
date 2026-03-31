import os

from gws_core import BaseTestCase, File, Folder, Settings, TaskRunner

from gws_omix import FastqFolder
from gws_omix.rna_seq.featurecounts.featurecounts import FeatureCounts
from gws_omix.rna_seq.mapping_genome.hisat2 import Hisat2Align
from gws_omix.rna_seq.mapping_genome.hisat2_index import Hisat2Index

TESTDATA = os.path.join(Settings.get_instance().get_variable("gws_omix", "testdata_dir"), "rna_seq")


class TestHisat2(BaseTestCase):

    def test_hisat2_index(self):
        """Hisat2Index builds a HISAT2 index from mini FASTA."""
        self.print("Test Hisat2Index: mini_genome.fa")

        output: Folder = TaskRunner(
            task_type=Hisat2Index,
            inputs={"genome_file": File(os.path.join(TESTDATA, "mini_genome.fa"))},
            params={"threads": 2},
        ).run()["output"]

        self.assertIsNotNone(output)
        ht2_files = [f for f in os.listdir(output.path) if f.startswith("genome_index")]
        self.print(f"HISAT2 index files: {ht2_files}")
        self.assertGreater(len(ht2_files), 0, msg="Expected .ht2 index files")

    def test_hisat2_align_se(self):
        """Chain Hisat2Index → Hisat2Align SE: produces BAM files."""
        self.print("Test Hisat2Align SE: build index then align mini_SE.fastq.gz")

        idx: Folder = TaskRunner(
            task_type=Hisat2Index,
            inputs={"genome_file": File(os.path.join(TESTDATA, "mini_genome.fa"))},
            params={"threads": 2},
        ).run()["output"]

        bam_folder: Folder = TaskRunner(
            task_type=Hisat2Align,
            inputs={
                "fastq_folder": FastqFolder(TESTDATA),
                "metadata": File(os.path.join(TESTDATA, "metadata_se.tsv")),
                "hisat2_index": idx,
            },
            params={"threads": 2},
        ).run()["output"]

        self.assertIsNotNone(bam_folder)
        bam_files = [f for f in os.listdir(bam_folder.path) if f.endswith(".bam")]
        self.print(f"BAM files: {bam_files}")
        self.assertGreater(len(bam_files), 0, msg="Expected BAM output from HISAT2")

    def test_featurecounts_chain(self):
        """Chain Hisat2Index → Hisat2Align (SE) → FeatureCounts: produces counts CSV."""
        self.print("Test FeatureCounts chain: index → align → count")

        idx: Folder = TaskRunner(
            task_type=Hisat2Index,
            inputs={"genome_file": File(os.path.join(TESTDATA, "mini_genome.fa"))},
            params={"threads": 2},
        ).run()["output"]

        bam_folder: Folder = TaskRunner(
            task_type=Hisat2Align,
            inputs={
                "fastq_folder": FastqFolder(TESTDATA),
                "metadata": File(os.path.join(TESTDATA, "metadata_se.tsv")),
                "hisat2_index": idx,
            },
            params={"threads": 2},
        ).run()["output"]

        counts_file: File = TaskRunner(
            task_type=FeatureCounts,
            inputs={
                "annotation_file": File(os.path.join(TESTDATA, "mini.gtf")),
                "bam_files": bam_folder,
            },
            params={"threads": 2, "sequencing_type": "Single-end", "strandedness": 0},
        ).run()["output"]

        self.assertIsNotNone(counts_file)
        self.assertTrue(os.path.exists(counts_file.path))
        import pandas as pd
        df = pd.read_csv(counts_file.path)
        self.print(f"Counts shape: {df.shape}")
        self.print(df.to_string())
        self.assertGreater(df.shape[0], 0, msg="Expected at least one gene in counts matrix")
