import os

from gws_core import BaseTestCase, File, Settings, TaskRunner

from gws_omix.phylogenetic_analysis.iqtree3.iqtree3 import IQTree3
from gws_omix.phylogenetic_analysis.phyTreeViz.phyTreeViz import PhyTreeViz

settings = Settings.get_instance()

# ---------------------------------------------------------------------------
# ADAPT THIS: aligned FASTA used to build the tree for visualization tests
# ---------------------------------------------------------------------------
ALIGNED_DNA_FASTA = "aligned_dna.fasta"


class TestPhyTreeViz(BaseTestCase):

    def _data_dir(self) -> str:
        return os.path.join(
            settings.get_variable("gws_omix", "testdata_dir"),
            "phylogenetic_analysis",
        )

    def _build_tree(self, bootstrap: int = 0) -> File:
        """Run IQTree3 to get a Newick tree file used as input for PhyTreeViz."""
        fa_path = os.path.join(self._data_dir(), ALIGNED_DNA_FASTA)
        runner = TaskRunner(
            task_type=IQTree3,
            inputs={"aligned_fasta": File(fa_path)},
            params={"prefix": "", "model": "MFP", "bootstrap": bootstrap, "threads": 1},
        )
        outputs = runner.run()
        return outputs["tree_file"]

    def _run_phytreeviz(self, tree_file: File, **extra_params):
        params = {
            "prefix": "",
            "fig_height": 3,
            "fig_width": 120,
            "leaf_label_size": 8,
            "align_leaf_label": True,
            "show_branch_length": False,
            "show_confidence": True,
            "dpi": 150,
        }
        params.update(extra_params)
        runner = TaskRunner(
            task_type=PhyTreeViz,
            inputs={"tree_file": tree_file},
            params=params,
        )
        return runner.run()

    # ------------------------------------------------------------------
    # Full chain: IQTree3 -> PhyTreeViz
    # ------------------------------------------------------------------

    def test_phytreeviz_chain_no_bootstrap(self):
        """Full chain IQTree3 (no bootstrap) -> PhyTreeViz -> PNG produced."""
        self.print("Test PhyTreeViz: chain IQTree3 (no bootstrap) -> PhyTreeViz")

        tree_file = self._build_tree(bootstrap=0)
        outputs = self._run_phytreeviz(tree_file)

        png = outputs["tree_figure"]
        self.assertTrue(os.path.isfile(png.path), "PNG figure not found")
        self.assertTrue(png.path.endswith(".png"))
        self.print(f"PNG: {png.path}")

    def test_phytreeviz_chain_with_bootstrap(self):
        """Full chain IQTree3 (bootstrap=1000) -> PhyTreeViz with confidence values shown."""
        self.print("Test PhyTreeViz: chain IQTree3 (bootstrap) -> PhyTreeViz + confidence")

        tree_file = self._build_tree(bootstrap=1000)
        outputs = self._run_phytreeviz(tree_file, show_confidence=True)

        png = outputs["tree_figure"]
        self.assertTrue(os.path.isfile(png.path), "PNG figure not found")

    # ------------------------------------------------------------------
    # Appearance options
    # ------------------------------------------------------------------

    def test_phytreeviz_show_branch_length(self):
        """show_branch_length=True -> PNG produced without error."""
        self.print("Test PhyTreeViz: show_branch_length=True")

        tree_file = self._build_tree()
        outputs = self._run_phytreeviz(tree_file, show_branch_length=True)

        self.assertTrue(os.path.isfile(outputs["tree_figure"].path))

    def test_phytreeviz_no_align_label(self):
        """align_leaf_label=False -> PNG produced without error."""
        self.print("Test PhyTreeViz: align_leaf_label=False")

        tree_file = self._build_tree()
        outputs = self._run_phytreeviz(tree_file, align_leaf_label=False)

        self.assertTrue(os.path.isfile(outputs["tree_figure"].path))

    def test_phytreeviz_custom_prefix(self):
        """Custom prefix -> output PNG filename contains the prefix."""
        self.print("Test PhyTreeViz: custom prefix")

        tree_file = self._build_tree()
        outputs = self._run_phytreeviz(tree_file, prefix="myfig")

        png = outputs["tree_figure"]
        self.assertTrue(os.path.isfile(png.path))
        self.assertIn("myfig", os.path.basename(png.path))
