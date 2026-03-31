import os

from gws_core import BaseTestCase, File, Settings, TaskRunner

from gws_omix.go_enrichment.gsea_analysis.gaf_to_gmt import GafToGmt
from gws_omix.go_enrichment.gsea_analysis.gsea import GSEAEnrichmentGMT

settings = Settings.get_instance()

# ---------------------------------------------------------------------------
# ADAPT THESE: file names and column names for your real data
# ---------------------------------------------------------------------------

# 1) GAF file (in testdata/go_enrichment/gaf_to_gmt/)
GAF_FILENAME = "sgd_gsea.gaf.gz"       # SGD yeast GAF (gzip)
GAF_PREFIX = "sgd"                     # prefix used to build GMT names

# The key from gmt_withIEA_by_go ResourceSet to use as GSEA gene sets.
# Generated keys are: "BP_id", "BP_symbol", "MF_id", "MF_symbol", "CC_id", "CC_symbol"
GMT_KEY = "BP_symbol"                  # SGD ORF names (YJR150C style) are in "symbol" column

# 2) DE results CSV (in testdata/go_enrichment/gsea/)
DE_FILENAME = "GSEA_DE_FULL_ALL_vs_CTRL.csv"

# Column names inside the DE CSV
GENES_COL = "gene_id"                  # column containing SGD gene IDs
STAT_COL = "stat"                      # column containing the signed ranking statistic


class TestGSEAEnrichmentGMT(BaseTestCase):

    def _gaf_dir(self) -> str:
        return os.path.join(
            settings.get_variable("gws_omix", "testdata_dir"),
            "go_enrichment", "gaf_to_gmt",
        )

    def _gsea_dir(self) -> str:
        return os.path.join(
            settings.get_variable("gws_omix", "testdata_dir"),
            "go_enrichment", "gsea",
        )

    def _build_gmt_from_gaf(self) -> File:
        """Run GafToGmt and return the GMT File for GMT_KEY (e.g. BP_symbol)."""
        gaf_path = os.path.join(self._gaf_dir(), GAF_FILENAME)
        runner = TaskRunner(
            task_type=GafToGmt,
            inputs={"gaf_file": File(gaf_path)},
            params={"prefix": GAF_PREFIX},
        )
        outputs = runner.run()
        gmt_resources = outputs["gmt_withIEA_by_go"].get_resources()
        self.assertIn(
            GMT_KEY, gmt_resources,
            f"GMT key '{GMT_KEY}' not found in gmt_withIEA_by_go. "
            f"Available keys: {list(gmt_resources.keys())}",
        )
        return gmt_resources[GMT_KEY]

    def _run_gsea(self, de_file: File, gmt_file: File, **extra_params):
        params = {
            "genes_colname": GENES_COL,
            "stat_colname": STAT_COL,
            "min_size": 5,
            "max_size": 500,
            "nperm": 1000,
            "processes": 1,
            "n_enrichplots": 3,
            "grid_genes_per_page": 25,
            "grid_sets_per_page": 10,
        }
        params.update(extra_params)
        runner = TaskRunner(
            task_type=GSEAEnrichmentGMT,
            inputs={"de_table_file": de_file, "gmt_file": gmt_file},
            params=params,
        )
        return runner.run()

    # ------------------------------------------------------------------
    # Tests
    # ------------------------------------------------------------------

    def test_gsea_from_gaf_gmt_chain(self):
        """Full chain: GAF -> GMT (GafToGmt) -> GSEA -> results table produced."""
        self.print("Test GSEAEnrichmentGMT: full chain GAF -> GMT -> GSEA")
        de_path = os.path.join(self._gsea_dir(), DE_FILENAME)

        gmt_file = self._build_gmt_from_gaf()
        outputs = self._run_gsea(File(de_path), gmt_file)

        results = outputs["results"]
        self.assertIsNotNone(results)
        rs_keys = results.get_resources()
        self.print(f"Results keys: {list(rs_keys.keys())}")
        self.assertIn("GSEA_results.csv", rs_keys)

    def test_gsea_enrichplots_output(self):
        """GSEA run -> enrichplots ResourceSet is produced."""
        self.print("Test GSEAEnrichmentGMT: enrichplots output")
        de_path = os.path.join(self._gsea_dir(), DE_FILENAME)

        gmt_file = self._build_gmt_from_gaf()
        outputs = self._run_gsea(File(de_path), gmt_file, n_enrichplots=3)

        enrichplots = outputs["enrichplots"]
        self.assertIsNotNone(enrichplots)
        self.print(f"Enrichplots keys: {list(enrichplots.get_resources().keys())}")

    def test_gsea_no_enrichplots_requested(self):
        """n_enrichplots=0 -> enrichplots ResourceSet is empty."""
        self.print("Test GSEAEnrichmentGMT: n_enrichplots=0")
        de_path = os.path.join(self._gsea_dir(), DE_FILENAME)

        gmt_file = self._build_gmt_from_gaf()
        outputs = self._run_gsea(File(de_path), gmt_file, n_enrichplots=0)

        enrichplots = outputs["enrichplots"]
        self.assertEqual(len(enrichplots.get_resources()), 0)
