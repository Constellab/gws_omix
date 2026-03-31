import os

from gws_core import BaseTestCase, File, Settings, TaskRunner

from gws_omix.go_enrichment.ora_analysis.ora_analysis import ORAEnrichmentTask

settings = Settings.get_instance()

# ---------------------------------------------------------------------------
# ADAPT THIS: file name of your DE results in testdata/go_enrichment/ora/
# ---------------------------------------------------------------------------
DE_FILENAME = "ORA_DE_FILT_ALL_vs_CTRL.csv"  # yeast filtered DE results
GENES_COL = "gene_id"                          # SGD gene ID column
PADJ_COL = "padj"                              # adjusted p-value column (used by task)
LOG2FC_COL = "log2FoldChange"                  # log2FC column name (informational)
ORGANISM = "Saccharomyces cerevisiae S288C"    # scientific name for g:Profiler


class TestORAEnrichmentTask(BaseTestCase):

    def _data_dir(self) -> str:
        return os.path.join(
            settings.get_variable("gws_omix", "testdata_dir"),
            "go_enrichment", "ora",
        )

    def _run_ora(self, de_file: File, correction_method: str = "gSCS", **extra_params):
        params = {
            "organism_name": ORGANISM,
            "genes_colname": GENES_COL,
            "padj_threshold": 0.05,
            "abs_log2fc": 0.0,
            "topn_plot": 10.0,
            "sources_list": "GO:BP,GO:MF,GO:CC",
            "correction_method": correction_method,
            "grid_genes_per_page": 25.0,
            "grid_terms_per_page": 20.0,
        }
        params.update(extra_params)
        runner = TaskRunner(
            task_type=ORAEnrichmentTask,
            inputs={"de_table_file": de_file},
            params=params,
        )
        return runner.run()

    # ------------------------------------------------------------------
    # Correction methods
    # ------------------------------------------------------------------

    def test_ora_gscs_correction(self):
        """ORA with gSCS correction (default g:Profiler method) -> tables, barplots, grids produced."""
        self.print("Test ORAEnrichmentTask: S. cerevisiae / gSCS correction")
        de_path = os.path.join(self._data_dir(), DE_FILENAME)

        outputs = self._run_ora(File(de_path), correction_method="gSCS")

        tables = outputs["tables"]
        barplots = outputs["barplots"]
        grids = outputs["term_gene_grid"]
        self.assertIsNotNone(tables)
        self.assertIsNotNone(barplots)
        self.assertIsNotNone(grids)
        self.print(f"Tables: {list(tables.get_resources().keys())}")
        self.print(f"Barplots: {list(barplots.get_resources().keys())}")
        self.print(f"Grids: {list(grids.get_resources().keys())}")

    def test_ora_fdr_correction(self):
        """ORA with Benjamini–Hochberg (fdr) correction -> outputs produced."""
        self.print("Test ORAEnrichmentTask: S. cerevisiae / fdr correction")
        de_path = os.path.join(self._data_dir(), DE_FILENAME)

        outputs = self._run_ora(File(de_path), correction_method="fdr")

        self.assertIsNotNone(outputs["tables"])
        self.assertIsNotNone(outputs["barplots"])
        self.assertIsNotNone(outputs["term_gene_grid"])

    def test_ora_bonferroni_correction(self):
        """ORA with Bonferroni correction -> outputs produced."""
        self.print("Test ORAEnrichmentTask: S. cerevisiae / bonferroni correction")
        de_path = os.path.join(self._data_dir(), DE_FILENAME)

        outputs = self._run_ora(File(de_path), correction_method="bonferroni")

        self.assertIsNotNone(outputs["tables"])
        self.assertIsNotNone(outputs["barplots"])
        self.assertIsNotNone(outputs["term_gene_grid"])

    # ------------------------------------------------------------------
    # log2FC threshold filtering
    # ------------------------------------------------------------------

    def test_ora_with_log2fc_threshold(self):
        """abs_log2fc=1.0 -> only genes with |log2FC| >= 1 are sent to g:Profiler."""
        self.print("Test ORAEnrichmentTask: abs_log2fc threshold = 1.0")
        de_path = os.path.join(self._data_dir(), DE_FILENAME)

        outputs = self._run_ora(
            File(de_path),
            correction_method="gSCS",
            abs_log2fc=1.0,
        )

        self.assertIsNotNone(outputs["tables"])
        self.assertIsNotNone(outputs["term_gene_grid"])

    # ------------------------------------------------------------------
    # Sources selection
    # ------------------------------------------------------------------

    def test_ora_kegg_only(self):
        """sources_list='KEGG' -> ORA restricted to KEGG pathways."""
        self.print("Test ORAEnrichmentTask: KEGG source only")
        de_path = os.path.join(self._data_dir(), DE_FILENAME)

        outputs = self._run_ora(
            File(de_path),
            correction_method="gSCS",
            sources_list="KEGG",
        )

        self.assertIsNotNone(outputs["tables"])
        self.assertIsNotNone(outputs["term_gene_grid"])
