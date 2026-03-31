import os

from gws_core import BaseTestCase, File, PlotlyResource, ResourceSet, Settings, TaskRunner

from gws_omix.rna_seq.multipydesq2 import Pydesq2Multi

TESTDATA = os.path.join(Settings.get_instance().get_variable("gws_omix", "testdata_dir"), "rna_seq")


class TestPydeseq2Multi(BaseTestCase):

    def test_deseq2_treat_vs_ctrl(self):
        """pyDESeq2 multi-contrast on 20-gene synthetic data (3 CTRL vs 3 TREAT)."""
        self.print("Test pyDESeq2: mini_counts.csv — TREAT vs CTRL")

        outputs = TaskRunner(
            task_type=Pydesq2Multi,
            inputs={
                "count_table_file": File(os.path.join(TESTDATA, "mini_counts.csv")),
                "metadata_file":    File(os.path.join(TESTDATA, "mini_metadata_deseq.tsv")),
                "gtf_file":         File(os.path.join(TESTDATA, "mini.gtf")),
            },
            params={
                "samples_colname":        "Sample",
                "genes_colname":          "gene_id",
                "condition_column":       "Condition",
                "control_condition":      "CTRL",
                "timepoint_column":       "",
                "group_column":           "",
                "extra_covariates":       "",
                "padj_value":             0.05,
                "log2FoldChange_value":   0.5,
                "annotation_columns":     "",
            },
        ).run()

        # --- all 5 output keys must be present ---
        for key in ("tables_full", "tables_filtered", "pca_plots", "heatmap_plots", "volcano_plots"):
            self.assertIn(key, outputs, msg=f"Missing output key: {key}")
            self.assertIsNotNone(outputs[key])

        # --- tables_full: at least 1 contrast ---
        tables_full: ResourceSet = outputs["tables_full"]
        resources = tables_full.get_resources()
        self.assertGreater(len(resources), 0, msg="Expected at least 1 DE results table")

        first_name, first_table = next(iter(resources.items()))
        df = first_table.get_data()
        self.print(f"Contrast '{first_name}': {df.shape[0]} genes × {df.shape[1]} cols")
        self.print(df.head(5).to_string())
        self.assertGreater(df.shape[0], 0, msg="DE results table should not be empty")

        # --- pca_plots: at least 1 PCA ---
        pca: ResourceSet = outputs["pca_plots"]
        self.assertGreater(len(pca.get_resources()), 0, msg="Expected PCA plot")
        self.print(f"PCA plots: {len(pca.get_resources())}")

        # --- volcano: at least 1 volcano per contrast ---
        volc: ResourceSet = outputs["volcano_plots"]
        self.assertGreater(len(volc.get_resources()), 0, msg="Expected volcano plot")
        self.print(f"Volcano plots: {len(volc.get_resources())}")

    def test_deseq2_with_timepoint_and_group(self):
        """pyDESeq2 multi-contrast with timepoint and group covariates."""
        self.print("Test pyDESeq2: mini_counts.csv — TREAT vs CTRL + timepoint + group")

        outputs = TaskRunner(
            task_type=Pydesq2Multi,
            inputs={
                "count_table_file": File(os.path.join(TESTDATA, "mini_counts.csv")),
                "metadata_file":    File(os.path.join(TESTDATA, "mini_metadata_deseq.tsv")),
                "gtf_file":         File(os.path.join(TESTDATA, "mini.gtf")),
            },
            params={
                "samples_colname":        "Sample",
                "genes_colname":          "gene_id",
                "condition_column":       "Condition",
                "control_condition":      "CTRL",
                "timepoint_column":       "timepoint",
                "group_column":           "group",
                "extra_covariates":       "",
                "padj_value":             0.05,
                "log2FoldChange_value":   0.5,
                "annotation_columns":     "",
            },
        ).run()

        for key in ("tables_full", "tables_filtered", "pca_plots", "heatmap_plots", "volcano_plots"):
            self.assertIn(key, outputs, msg=f"Missing output key: {key}")
            self.assertIsNotNone(outputs[key])

        tables_full: ResourceSet = outputs["tables_full"]
        resources = tables_full.get_resources()
        self.assertGreater(len(resources), 0, msg="Expected at least 1 DE results table")
        self.print(f"Contrasts with timepoint+group: {list(resources.keys())}")

        pca: ResourceSet = outputs["pca_plots"]
        self.assertGreater(len(pca.get_resources()), 0, msg="Expected PCA plot")
