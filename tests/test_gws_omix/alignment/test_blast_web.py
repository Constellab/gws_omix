import os

from gws_core import BaseTestCase, File, ResourceSet, Settings, Table, TaskRunner

from gws_omix.aligment.blast_web import BlastWebRunner


def _count_sequences(fasta_path: str) -> int:
    with open(fasta_path) as f:
        return sum(1 for line in f if line.startswith(">"))

settings = Settings.get_instance()

# ---------------------------------------------------------------------------
# Testdata
# ---------------------------------------------------------------------------
TESTDATA_DIR = settings.get_variable("gws_omix", "testdata_dir")
PROTEIN_FASTA = os.path.join(TESTDATA_DIR, "alignment", "small_proteins.fasta")
NUCLEOTIDE_FASTA = os.path.join(TESTDATA_DIR, "alignment", "small_nucleotides.fasta")


class TestBlastWebRunner(BaseTestCase):

    def _run(self, fasta_path: str, blast_program: str, sequence_type: str,
             evalue: float = 1e-3, max_target_seqs: int = 5) -> dict:
        runner = TaskRunner(
            task_type=BlastWebRunner,
            inputs={"input_fasta": File(fasta_path)},
            params={
                "blast_program": blast_program,
                "sequence_type": sequence_type,
                "evalue": evalue,
                "max_target_seqs": max_target_seqs,
            },
        )
        return runner.run()

    # ------------------------------------------------------------------
    # blastn — nucleotide query vs NCBI nt (remote)
    # ------------------------------------------------------------------

    def test_blastn(self):
        """Nucleotide FASTA → remote blastn → ResourceSet of result tables."""
        input_seqs = _count_sequences(NUCLEOTIDE_FASTA)
        self.print(f"Test BlastWebRunner: blastn (remote NCBI)  |  Input sequences: {input_seqs}")

        outputs = self._run(NUCLEOTIDE_FASTA, blast_program="blastn", sequence_type="nucl")

        results: ResourceSet = outputs["blast_results"]
        self.assertIsNotNone(results)
        self.print(f"Number of result tables: {len(results)}  (one per input sequence)")
        for name, resource in results.get_resources().items():
            df = resource.get_data()
            self.print(f"\n--- Table: {name}  ({df.shape[0]} rows x {df.shape[1]} cols) ---")
            self.print(df.to_string())

    # ------------------------------------------------------------------
    # blastp — protein query vs NCBI nr (remote)
    # ------------------------------------------------------------------

    def test_blastp(self):
        """Protein FASTA → remote blastp → ResourceSet of result tables."""
        input_seqs = _count_sequences(PROTEIN_FASTA)
        self.print(f"Test BlastWebRunner: blastp (remote NCBI)  |  Input sequences: {input_seqs}")

        outputs = self._run(PROTEIN_FASTA, blast_program="blastp", sequence_type="prot")

        results: ResourceSet = outputs["blast_results"]
        self.assertIsNotNone(results)
        self.print(f"Number of result tables: {len(results)}  (one per input sequence)")
        for name, resource in results.get_resources().items():
            df = resource.get_data()
            self.print(f"\n--- Table: {name}  ({df.shape[0]} rows x {df.shape[1]} cols) ---")
            self.print(df.to_string())
