import os

import pandas as pd
from gws_core import BaseTestCase, File, Settings, TaskRunner

from gws_omix.file.genes_id_conversion.genes_id_conversion import IDConvertTask

settings = Settings.get_instance()

# ---------------------------------------------------------------------------
# Input file: reuse the yeast ORA DE results (contains gene_id column)
# Organism: Saccharomyces cerevisiae S288C
# ---------------------------------------------------------------------------
DE_FILE = os.path.join(
    settings.get_variable("gws_omix", "testdata_dir"),
    "go_enrichment", "ora", "ORA_DE_FILT_ALL_vs_CTRL.csv",
)
ORGANISM = "Saccharomyces cerevisiae S288C"
ID_COLUMN = "gene_id"   # SGD ORF names (YJR150C style)


class TestIDConvert(BaseTestCase):

    def _run(self, target_namespace: str, **extra_params):
        params = {
            "organism_name": ORGANISM,
            "id_column": ID_COLUMN,
            "target_namespace": target_namespace,
            "numeric_namespace": "auto",
        }
        params.update(extra_params)
        runner = TaskRunner(
            task_type=IDConvertTask,
            inputs={"table_file": File(DE_FILE)},
            params=params,
        )
        return runner.run()

    # ------------------------------------------------------------------
    # Target namespace: ENSG (Ensembl gene IDs)
    # ------------------------------------------------------------------

    def test_convert_to_ensg(self):
        """SGD ORF names -> Ensembl gene IDs (ENSG) -> converted_table produced."""
        self.print("Test IDConvert: gene_id -> ENSG")

        input_rows = pd.read_csv(DE_FILE).shape[0]
        self.print(f"Input rows:     {input_rows}")

        outputs = self._run("ENSG")

        converted = outputs["converted_table"]
        annotated = outputs["annotated_file"]
        self.assertIsNotNone(converted)
        self.assertIsNotNone(annotated)

        converted_rows = converted.get_data().shape[0]
        annotated_rows = annotated.get_data().shape[0]
        self.print(f"Converted rows: {converted_rows}  (input: {input_rows})")
        self.print(f"Annotated rows: {annotated_rows}  (input: {input_rows})")

    # ------------------------------------------------------------------
    # Target namespace: UNIPROTSWISSPROT
    # ------------------------------------------------------------------

    def test_convert_to_uniprot(self):
        """SGD ORF names -> UniProtKB/Swiss-Prot accessions."""
        self.print("Test IDConvert: gene_id -> UNIPROTSWISSPROT")

        outputs = self._run("UNIPROTSWISSPROT")

        self.assertIsNotNone(outputs["converted_table"])
        self.assertIsNotNone(outputs["annotated_file"])

    # ------------------------------------------------------------------
    # Target namespace: ENTREZGENE
    # ------------------------------------------------------------------

    def test_convert_to_entrezgene(self):
        """SGD ORF names -> NCBI Entrez Gene IDs."""
        self.print("Test IDConvert: gene_id -> ENTREZGENE")

        outputs = self._run("ENTREZGENE", numeric_namespace="ENTREZGENE_ACC")

        self.assertIsNotNone(outputs["converted_table"])
        self.assertIsNotNone(outputs["annotated_file"])

    # ------------------------------------------------------------------
    # annotated_file always contains the original input rows
    # ------------------------------------------------------------------

    def test_annotated_file_contains_original_column(self):
        """annotated_file must keep the original gene_id column."""
        self.print("Test IDConvert: annotated_file keeps original gene_id column")

        outputs = self._run("ENSG")

        annotated = outputs["annotated_file"]
        self.assertIsNotNone(annotated)
        col_names = list(annotated.get_data().columns)
        self.print(f"Annotated columns: {col_names}")
        # The original gene_id column must survive in the annotated table
        self.assertIn(ID_COLUMN, col_names)
