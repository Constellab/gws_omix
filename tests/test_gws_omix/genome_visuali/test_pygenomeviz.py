import os

from gws_core import BaseTestCase, File, Folder, Settings, TaskRunner

from gws_omix.genome_visuali.pyGenomeViz import PyGenomeVizCompare

settings = Settings.get_instance()


class TestPyGenomeVizCompare(BaseTestCase):

    def _data_dir(self) -> str:
        return os.path.join(
            settings.get_variable("gws_omix", "testdata_dir"),
            "genome_visuali",
        )

    def _run(self, genome_input, comparison_mode: str, prefix: str = "pygenomeviz"):
        runner = TaskRunner(
            task_type=PyGenomeVizCompare,
            inputs={"genome_input": genome_input},
            params={
                "prefix": prefix,
                "comparison_mode": comparison_mode,
            },
        )
        return runner.run()

    # ------------------------------------------------------------------
    # visualization_only mode
    # ------------------------------------------------------------------

    def test_visualization_only_single_genbank(self):
        """visualization_only with a single GenBank file -> HTML report."""
        self.print("Test PyGenomeVizCompare: visualization_only / single GenBank")
        gb_path = os.path.join(self._data_dir(), "ecoli_NC002483.gb")

        outputs = self._run(File(gb_path), comparison_mode="visualization_only")

        html = outputs["report_html"]
        self.assertTrue(os.path.isfile(html.path))
        self.assertTrue(html.path.endswith(".html"))

    def test_visualization_only_single_gff3(self):
        """visualization_only with a single GFF3 file -> HTML report."""
        self.print("Test PyGenomeVizCompare: visualization_only / single GFF3")
        gff_path = os.path.join(self._data_dir(), "ecoli_NC002483.gff3")

        outputs = self._run(File(gff_path), comparison_mode="visualization_only")

        html = outputs["report_html"]
        self.assertTrue(os.path.isfile(html.path))
        self.assertTrue(html.path.endswith(".html"))

    def test_visualization_only_folder_multiple_genbank(self):
        """visualization_only with a folder of 4 GenBank files -> HTML report."""
        self.print("Test PyGenomeVizCompare: visualization_only / folder 4 GenBank")
        folder_path = os.path.join(self._data_dir(), "GenomeVizCompare")

        outputs = self._run(Folder(folder_path), comparison_mode="visualization_only")

        html = outputs["report_html"]
        self.assertTrue(os.path.isfile(html.path))
        self.assertTrue(html.path.endswith(".html"))

    # ------------------------------------------------------------------
    # genome_comparison mode  (requires nucmer / pgv-mummer in PATH)
    # ------------------------------------------------------------------

    def test_genome_comparison_multi_genbank(self):
        """genome_comparison with 4 GenBank files -> HTML report with NUCmer links."""
        self.print("Test PyGenomeVizCompare: genome_comparison / 4 GenBank files")
        folder_path = os.path.join(self._data_dir(), "GenomeVizCompare")

        outputs = self._run(Folder(folder_path), comparison_mode="genome_comparison")

        html = outputs["report_html"]
        self.assertTrue(os.path.isfile(html.path))
        self.assertTrue(html.path.endswith(".html"))

    # ------------------------------------------------------------------
    # cds_protein_homology mode  (requires mmseqs / pgv-mmseqs in PATH)
    # ------------------------------------------------------------------

    def test_cds_protein_homology_multi_genbank(self):
        """cds_protein_homology with 4 GenBank files -> HTML report with MMseqs2 links."""
        self.print("Test PyGenomeVizCompare: cds_protein_homology / 4 GenBank files")
        folder_path = os.path.join(self._data_dir(), "GenomeVizCompare")

        outputs = self._run(
            Folder(folder_path),
            comparison_mode="cds_protein_homology",
        )

        html = outputs["report_html"]
        self.assertTrue(os.path.isfile(html.path))
        self.assertTrue(html.path.endswith(".html"))
