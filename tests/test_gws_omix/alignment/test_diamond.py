import os

from gws_core import BaseTestCase, File, Settings, Table, TaskRunner


def _count_sequences(fasta_path: str) -> int:
    """Count the number of sequences in a FASTA file."""
    with open(fasta_path) as f:
        return sum(1 for line in f if line.startswith(">"))

from gws_omix.aligment.blast_diamond import Diamond
from gws_omix.aligment.diamond_ECnumber import ECnumber

settings = Settings.get_instance()

# ---------------------------------------------------------------------------
# Testdata
# ---------------------------------------------------------------------------
TESTDATA_DIR = settings.get_variable("gws_omix", "testdata_dir")
PROTEIN_FASTA = os.path.join(TESTDATA_DIR, "alignment", "small_proteins.fasta")
NUCLEOTIDE_FASTA = os.path.join(TESTDATA_DIR, "alignment", "small_nucleotides.fasta")

# Diamond output columns (set by the task)
EXPECTED_COLUMNS = [
    "Subject ID", "Percentage of identical matches", "Alignment length",
    "Number of mismatches", "Number of gap openings",
    "Start of alignment in query", "End of alignment in query",
    "Start of alignment in subject", "End of alignment in subject",
    "Expected value", "Bit score",
]


class TestDiamond(BaseTestCase):

    def _run(self, fasta_path: str, input_type: str = "prot", evalue: float = 1e-5,
             threads: int = 1, query_cover: int = 70) -> dict:
        runner = TaskRunner(
            task_type=Diamond,
            inputs={"input_path": File(fasta_path)},
            params={
                "input_type_value": input_type,
                "evalue_value": evalue,
                "num_threads": threads,
                "query_cover_value": query_cover,
            },
        )
        return runner.run()

    # ------------------------------------------------------------------
    # blastp mode (protein query vs SwissProt protein DB)
    # ------------------------------------------------------------------

    def test_blastp(self):
        """Protein FASTA → Diamond blastp → Table with 11 expected columns."""
        self.print("Test Diamond: blastp (prot vs SwissProt)")

        input_seqs = _count_sequences(PROTEIN_FASTA)
        self.print(f"Input sequences:  {input_seqs}")

        outputs = self._run(PROTEIN_FASTA, input_type="prot")

        result: Table = outputs["output_path"]
        self.assertIsNotNone(result)
        output_rows = result.get_data().shape[0]
        cols = list(result.get_data().columns)
        self.print(f"Output columns:   {cols}")
        self.print(f"Output rows:      {output_rows}  (input sequences: {input_seqs}, avg hits/seq: {output_rows / input_seqs:.1f})")
        self.print("\n" + result.get_data().to_string())
        for col in EXPECTED_COLUMNS:
            self.assertIn(col, cols, msg=f"Missing column: {col}")

    # ------------------------------------------------------------------
    # blastx mode (nucleotide query translated vs SwissProt protein DB)
    # ------------------------------------------------------------------

    def test_blastx(self):
        """Nucleotide FASTA → Diamond blastx → Table with 11 expected columns."""
        self.print("Test Diamond: blastx (nuc translated vs SwissProt)")

        input_seqs = _count_sequences(NUCLEOTIDE_FASTA)
        self.print(f"Input sequences:  {input_seqs}")

        outputs = self._run(NUCLEOTIDE_FASTA, input_type="nuc")

        result: Table = outputs["output_path"]
        self.assertIsNotNone(result)
        output_rows = result.get_data().shape[0]
        self.print(f"Output rows:      {output_rows}  (input sequences: {input_seqs}, avg hits/seq: {output_rows / input_seqs:.1f})")
        self.print("\n" + result.get_data().to_string())
        for col in EXPECTED_COLUMNS:
            self.assertIn(col, list(result.get_data().columns), msg=f"Missing column: {col}")

    # ------------------------------------------------------------------
    # Chain: Diamond (blastp) → ECnumber
    # ------------------------------------------------------------------

    def test_blastp_to_ecnumber(self):
        """Diamond blastp → ECnumber → output Table produced (may be empty if no fungi EC hit)."""
        self.print("Test Diamond → ECnumber chain")

        input_seqs = _count_sequences(PROTEIN_FASTA)
        self.print(f"Input sequences:  {input_seqs}")

        diamond_outputs = self._run(PROTEIN_FASTA, input_type="prot")
        diamond_table: Table = diamond_outputs["output_path"]
        diamond_rows = diamond_table.get_data().shape[0]
        self.print(f"Diamond rows:     {diamond_rows}  (input sequences: {input_seqs}, avg hits/seq: {diamond_rows / input_seqs:.1f})")

        ec_runner = TaskRunner(
            task_type=ECnumber,
            inputs={"input_table": diamond_table},
            params={},
        )
        ec_outputs = ec_runner.run()

        ec_table: Table = ec_outputs["output_table"]
        self.assertIsNotNone(ec_table)
        self.print(f"ECnumber output rows: {ec_table.get_data().shape[0]}")
        self.print(f"ECnumber output cols: {list(ec_table.get_data().columns)}")
        self.print("\n" + ec_table.get_data().to_string())
