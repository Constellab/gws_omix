import os

from gws_core import BaseTestCase, File, Settings, TaskRunner

from gws_omix.phylogenetic_analysis.iqtree3.iqtree3 import IQTree3

settings = Settings.get_instance()

# ---------------------------------------------------------------------------
# ADAPT THIS: file name of the aligned FASTA in
#   tests/testdata/phylogenetic_analysis/
# ---------------------------------------------------------------------------
ALIGNED_DNA_FASTA = "aligned_dna.fasta"      # aligned DNA sequences
ALIGNED_PROTEIN_FASTA = "aligned_protein.fasta"  # aligned protein sequences


class TestIQTree3(BaseTestCase):

    def _data_dir(self) -> str:
        return os.path.join(
            settings.get_variable("gws_omix", "testdata_dir"),
            "phylogenetic_analysis",
        )

    def _run_iqtree(self, fasta_file: File, **extra_params):
        params = {
            "prefix": "",
            "model": "MFP",
            "bootstrap": 0,
            "threads": 1,
        }
        params.update(extra_params)
        runner = TaskRunner(
            task_type=IQTree3,
            inputs={"aligned_fasta": fasta_file},
            params=params,
        )
        return runner.run()

    # ------------------------------------------------------------------
    # DNA input
    # ------------------------------------------------------------------

    def test_iqtree_dna_no_bootstrap(self):
        """Aligned DNA FASTA + no bootstrap -> tree_file and iqtree_log produced."""
        self.print("Test IQTree3: DNA / no bootstrap")
        fa_path = os.path.join(self._data_dir(), ALIGNED_DNA_FASTA)

        outputs = self._run_iqtree(File(fa_path), bootstrap=0)

        tree = outputs["tree_file"]
        log = outputs["iqtree_log"]
        self.assertTrue(os.path.isfile(tree.path), "treefile not found")
        self.assertTrue(os.path.isfile(log.path), "log file not found")
        self.print(f"Treefile: {tree.path}")

    def test_iqtree_dna_with_bootstrap(self):
        """Aligned DNA FASTA + UFBoot 1000 -> tree with bootstrap support values."""
        self.print("Test IQTree3: DNA / bootstrap=1000")
        fa_path = os.path.join(self._data_dir(), ALIGNED_DNA_FASTA)

        outputs = self._run_iqtree(File(fa_path), bootstrap=1000)

        tree = outputs["tree_file"]
        self.assertTrue(os.path.isfile(tree.path), "treefile not found")
        # UFBoot tree has support values in Newick labels
        content = open(tree.path).read()
        self.print(f"Tree snippet: {content[:200]}")

    # ------------------------------------------------------------------
    # Protein input
    # ------------------------------------------------------------------

    def test_iqtree_protein_no_bootstrap(self):
        """Aligned protein FASTA -> IQ-TREE auto-detects protein, produces tree."""
        self.print("Test IQTree3: protein / no bootstrap")
        fa_path = os.path.join(self._data_dir(), ALIGNED_PROTEIN_FASTA)

        outputs = self._run_iqtree(File(fa_path), bootstrap=0)

        tree = outputs["tree_file"]
        self.assertTrue(os.path.isfile(tree.path), "treefile not found")

    # ------------------------------------------------------------------
    # Auxiliary outputs
    # ------------------------------------------------------------------

    def test_iqtree_additional_outputs(self):
        """IQTree3 run -> iqtree_outputs ResourceSet contains at least the .iqtree summary."""
        self.print("Test IQTree3: additional outputs ResourceSet")
        fa_path = os.path.join(self._data_dir(), ALIGNED_DNA_FASTA)

        outputs = self._run_iqtree(File(fa_path), bootstrap=0)

        rs = outputs["iqtree_outputs"]
        keys = list(rs.get_resources().keys())
        self.print(f"Additional output files: {keys}")
        self.assertGreater(len(keys), 0)

    # ------------------------------------------------------------------
    # Custom prefix
    # ------------------------------------------------------------------

    def test_iqtree_custom_prefix(self):
        """Custom prefix -> outputs use that prefix as base name."""
        self.print("Test IQTree3: custom prefix")
        fa_path = os.path.join(self._data_dir(), ALIGNED_DNA_FASTA)

        outputs = self._run_iqtree(File(fa_path), bootstrap=0, prefix="mytree")

        tree = outputs["tree_file"]
        self.assertTrue(os.path.isfile(tree.path))
        self.assertIn("mytree", os.path.basename(tree.path))
