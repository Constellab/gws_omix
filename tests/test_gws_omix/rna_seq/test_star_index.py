import os

from gws_core import BaseTestCase, File, Folder, Settings, TaskRunner

from gws_omix.rna_seq.star_index.star_index import StarIndex

TESTDATA = os.path.join(Settings.get_instance().get_variable("gws_omix", "testdata_dir"), "rna_seq")


class TestStarIndex(BaseTestCase):

    def test_star_index_builds(self):
        """StarIndex builds a STAR genome index from mini FASTA + GTF."""
        self.print("Test StarIndex: mini_genome.fa + mini.gtf")

        output: Folder = TaskRunner(
            task_type=StarIndex,
            inputs={
                "genome_file": File(os.path.join(TESTDATA, "mini_genome.fa")),
                "annotation_file": File(os.path.join(TESTDATA, "mini.gtf")),
            },
            params={"threads": 2},
        ).run()["output"]

        self.assertIsNotNone(output)
        index_files = os.listdir(output.path)
        self.print(f"STAR index files: {index_files}")
        self.assertGreater(len(index_files), 0, msg="Expected STAR index files in output folder")
