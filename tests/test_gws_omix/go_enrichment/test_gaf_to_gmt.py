import os

from gws_core import BaseTestCase, File, Settings, TaskRunner

from gws_omix.go_enrichment.gsea_analysis.gaf_to_gmt import GafToGmt

settings = Settings.get_instance()

# ---------------------------------------------------------------------------
# ADAPT THIS: name of the GAF (or GAF.GZ) file to put in
#   tests/testdata/go_enrichment/gaf_to_gmt/
# ---------------------------------------------------------------------------
GAF_FILENAME = "sgd_gsea.gaf.gz"   # SGD yeast GAF (gzip)
GAF_PREFIX = "sgd"                 # prefix used for GMT output file names


class TestGafToGmt(BaseTestCase):

    def _data_dir(self) -> str:
        return os.path.join(
            settings.get_variable("gws_omix", "testdata_dir"),
            "go_enrichment", "gaf_to_gmt",
        )

    def _run_gaf_to_gmt(self, gaf_file: File, prefix: str = GAF_PREFIX):
        runner = TaskRunner(
            task_type=GafToGmt,
            inputs={"gaf_file": gaf_file},
            params={"prefix": prefix},
        )
        return runner.run()

    def test_gaf_to_gmt_outputs_produced(self):
        """GAF file -> at least one GMT file in each output ResourceSet."""
        self.print("Test GafToGmt: all three output ResourceSets produced")
        gaf_path = os.path.join(self._data_dir(), GAF_FILENAME)

        outputs = self._run_gaf_to_gmt(File(gaf_path))

        gmt_all = outputs["gmt_files"].get_resources()
        gmt_with = outputs["gmt_withIEA_by_go"].get_resources()
        gmt_no = outputs["gmt_noIEA_by_go"].get_resources()

        self.print(f"gmt_files keys:   {list(gmt_all.keys())}")
        self.print(f"gmt_withIEA keys: {list(gmt_with.keys())}")
        self.print(f"gmt_noIEA keys:   {list(gmt_no.keys())}")

        # GO_ALL (4 files: withIEA/noIEA × id/symbol)
        self.assertGreater(len(gmt_all), 0)
        self.assertIn("withIEA_id", gmt_all)
        self.assertIn("noIEA_id", gmt_all)

    def test_gaf_to_gmt_with_iea_split_by_aspect(self):
        """GAF file -> withIEA split ResourceSet contains BP, MF and/or CC GMT files."""
        self.print("Test GafToGmt: withIEA split by GO aspect")
        gaf_path = os.path.join(self._data_dir(), GAF_FILENAME)

        outputs = self._run_gaf_to_gmt(File(gaf_path))

        gmt_with = outputs["gmt_withIEA_by_go"].get_resources()
        found_aspects = {k.split("_")[0] for k in gmt_with.keys()}
        self.print(f"Aspects found (withIEA): {found_aspects}")
        self.assertTrue(found_aspects & {"BP", "MF", "CC"})

    def test_gaf_to_gmt_no_iea_split_by_aspect(self):
        """GAF file -> noIEA split ResourceSet is created (IEA annotations excluded)."""
        self.print("Test GafToGmt: noIEA split by GO aspect")
        gaf_path = os.path.join(self._data_dir(), GAF_FILENAME)

        outputs = self._run_gaf_to_gmt(File(gaf_path))

        gmt_no = outputs["gmt_noIEA_by_go"].get_resources()
        self.print(f"noIEA keys: {list(gmt_no.keys())}")
        self.assertGreater(len(gmt_no), 0)

    def test_gaf_to_gmt_custom_prefix(self):
        """Custom prefix -> task completes without error."""
        self.print("Test GafToGmt: custom prefix")
        gaf_path = os.path.join(self._data_dir(), GAF_FILENAME)

        outputs = self._run_gaf_to_gmt(File(gaf_path), prefix="myorganism")

        gmt_all = outputs["gmt_files"].get_resources()
        self.assertGreater(len(gmt_all), 0)
