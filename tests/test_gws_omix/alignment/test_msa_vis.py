import os

from gws_core import BaseTestCase, File, ResourceSet, Settings, TaskRunner

from gws_omix.aligment.multiple_sequence_alignment_vis import MSAVis


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
# Already aligned — used for align_if_needed=False test
ALIGNED_PROTEIN_FASTA = os.path.join(TESTDATA_DIR, "phylogenetic_analysis", "aligned_protein.fasta")


class TestMSAVis(BaseTestCase):

    def _run(self, fasta_path: str, **extra_params) -> dict:
        params = {
            "prefix": "",
            "align_if_needed": True,
            "threads": 4,
        }
        params.update(extra_params)
        runner = TaskRunner(
            task_type=MSAVis,
            inputs={"sequences_fasta": File(fasta_path)},
            params=params,
        )
        return runner.run()

    # ------------------------------------------------------------------
    # Protein sequences — auto MAFFT alignment + visualization
    # ------------------------------------------------------------------

    def test_protein_default(self):
        """Protein FASTA → MSAVis → PNGs + aligned FASTA + log."""
        self.print("Test MSAVis: protein sequences, default params")

        outputs = self._run(PROTEIN_FASTA)

        msa_outputs: ResourceSet = outputs["msa_outputs"]
        aligned_fasta: File = outputs["aligned_fasta"]
        alignment_log: File = outputs["alignment_log"]

        self.assertIsNotNone(msa_outputs)
        self.assertIsNotNone(aligned_fasta)
        self.assertIsNotNone(alignment_log)
        input_seqs = _count_sequences(PROTEIN_FASTA)
        self.print(f"Input sequences:  {input_seqs}")
        self.print(f"PNG outputs:      {len(msa_outputs)}")
        self.print(f"PNG names:        {list(msa_outputs.get_resources().keys())}")
        self.print(f"Aligned FASTA:    {aligned_fasta.path}")
        self.print(f"Alignment log:    {alignment_log.path}")
        with open(aligned_fasta.path) as f:
            self.print("\n--- Aligned FASTA (first 20 lines) ---")
            self.print("".join(f.readlines()[:20]))

    # ------------------------------------------------------------------
    # Nucleotide sequences
    # ------------------------------------------------------------------

    def test_nucleotide_default(self):
        """Nucleotide FASTA → MSAVis → PNGs + aligned FASTA + log."""
        self.print("Test MSAVis: nucleotide sequences")

        outputs = self._run(NUCLEOTIDE_FASTA)

        self.assertIsNotNone(outputs["msa_outputs"])
        self.assertIsNotNone(outputs["aligned_fasta"])
        self.assertIsNotNone(outputs["alignment_log"])
        input_seqs = _count_sequences(NUCLEOTIDE_FASTA)
        self.print(f"Input sequences: {input_seqs}")
        self.print(f"PNG outputs:     {len(outputs['msa_outputs'])}")
        self.print(f"PNG names:       {list(outputs['msa_outputs'].get_resources().keys())}")

    # ------------------------------------------------------------------
    # Custom prefix
    # ------------------------------------------------------------------

    def test_custom_prefix(self):
        """MSAVis with custom prefix → PNG filename starts with the prefix."""
        self.print("Test MSAVis: custom prefix 'myalign'")

        outputs = self._run(PROTEIN_FASTA, prefix="myalign")

        msa_outputs: ResourceSet = outputs["msa_outputs"]
        self.assertIsNotNone(msa_outputs)
        png_names = list(msa_outputs.get_resources().keys())
        self.print(f"PNG names: {png_names}")
        for name in png_names:
            self.assertTrue(
                name.startswith("myalign"),
                msg=f"Expected PNG name to start with 'myalign', got: {name}"
            )

    # ------------------------------------------------------------------
    # align_if_needed=False vs True — même fichier pré-aligné, comparaison directe
    # ------------------------------------------------------------------

    def test_align_if_needed_false_vs_true(self):
        """Same pre-aligned FASTA run with align_if_needed=False then True → compare outputs."""
        input_seqs = _count_sequences(ALIGNED_PROTEIN_FASTA)
        self.print(f"Input sequences: {input_seqs}")

        # --- Case False: skip MAFFT, visualize directly ---
        self.print("\n=== align_if_needed=False (direct visualization, no MAFFT) ===")
        out_false = self._run(ALIGNED_PROTEIN_FASTA, align_if_needed=False, prefix="no_mafft")
        pngs_false = list(out_false["msa_outputs"].get_resources().keys())
        self.assertIsNotNone(out_false["msa_outputs"])
        self.assertIsNotNone(out_false["aligned_fasta"])
        self.print(f"PNG outputs: {len(out_false['msa_outputs'])}")
        self.print(f"PNG names:   {pngs_false}")

        # --- Case True: run MAFFT first even though already aligned ---
        self.print("\n=== align_if_needed=True (MAFFT runs first, then visualization) ===")
        out_true = self._run(ALIGNED_PROTEIN_FASTA, align_if_needed=True, prefix="with_mafft")
        pngs_true = list(out_true["msa_outputs"].get_resources().keys())
        self.assertIsNotNone(out_true["msa_outputs"])
        self.assertIsNotNone(out_true["aligned_fasta"])
        self.print(f"PNG outputs: {len(out_true['msa_outputs'])}")
        self.print(f"PNG names:   {pngs_true}")

        # Both should produce the same number of PNGs
        self.assertEqual(
            len(out_false["msa_outputs"]),
            len(out_true["msa_outputs"]),
            msg="Both modes should produce the same number of PNG outputs"
        )
        self.print(f"\nBoth modes produced {len(pngs_false)} PNG(s) ✓")
