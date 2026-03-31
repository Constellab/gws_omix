import os

from gws_core import BaseTestCase, File, Folder, Settings, TaskRunner

from gws_omix.genome_visuali.pyCirclize import PyCirclize

settings = Settings.get_instance()


class TestPyCirclize(BaseTestCase):

    def _data_dir(self) -> str:
        return os.path.join(
            settings.get_variable("gws_omix", "testdata_dir"),
            "genome_visuali",
        )

    def _run_pycirclize(self, genome_input):
        runner = TaskRunner(
            task_type=PyCirclize,
            inputs={"genome_input": genome_input},
            params={},
        )
        return runner.run()

    def test_single_genbank_file(self):
        """Single GenBank .gb file -> ResourceSet containing exactly 1 PNG."""
        self.print("Test PyCirclize: single GenBank file")
        gb_path = os.path.join(self._data_dir(), "ecoli_NC002483.gb")

        outputs = self._run_pycirclize(File(gb_path))

        plots = outputs["circular_plots"].get_resources()
        self.print(f"PNGs generated: {list(plots.keys())}")
        self.assertEqual(len(plots), 1)
        self.assertIn("ecoli_NC002483.png", plots)

    def test_single_gff3_file(self):
        """Single GFF3 .gff3 file -> ResourceSet containing exactly 1 PNG."""
        self.print("Test PyCirclize: single GFF3 file")
        gff_path = os.path.join(self._data_dir(), "ecoli_NC002483.gff3")

        outputs = self._run_pycirclize(File(gff_path))

        plots = outputs["circular_plots"].get_resources()
        self.print(f"PNGs generated: {list(plots.keys())}")
        self.assertEqual(len(plots), 1)
        self.assertIn("ecoli_NC002483.png", plots)

    def test_folder_multiple_genbank_files(self):
        """Folder with 4 GenBank files -> ResourceSet containing 4 PNGs (one per genome)."""
        self.print("Test PyCirclize: folder with multiple GenBank files")
        folder_path = os.path.join(self._data_dir(), "GenomeVizCompare")

        outputs = self._run_pycirclize(Folder(folder_path))

        plots = outputs["circular_plots"].get_resources()
        self.print(f"PNGs generated: {list(plots.keys())}")
        self.assertEqual(len(plots), 4)
        self.assertIn("NC_070914.png", plots)
        self.assertIn("NC_070915.png", plots)
        self.assertIn("NC_070916.png", plots)
        self.assertIn("NC_070918.png", plots)
