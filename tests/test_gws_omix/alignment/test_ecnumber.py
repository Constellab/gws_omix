import os

import pandas as pd
from gws_core import BaseTestCase, File, Settings, Table, TaskRunner

from gws_omix.aligment.diamond_ECnumber import ECnumber

# ---------------------------------------------------------------------------
# Known S. cerevisiae (fungi) UniProt IDs with EC numbers in the DB.
# Diamond outputs Subject IDs in "ACCESSION.VERSION" format (e.g. P00330.5).
# The task strips the version suffix, so we simulate that format here.
# ---------------------------------------------------------------------------
SUBJECT_IDS = [
    "P00330.5",   # ADH1_YEAST — alcohol dehydrogenase   → EC 1.1.1.1
    "P00359.3",   # TDH3_YEAST — GAPDH                   → EC 1.2.1.12
    "P00549.1",   # PYK1_YEAST — pyruvate kinase          → EC 2.7.1.40
    "FAKE999.1",  # intentionally unknown — no EC hit expected
]


def _make_diamond_table(subject_ids: list) -> Table:
    """Build a minimal Diamond-like Table with a 'Subject ID' column."""
    df = pd.DataFrame({
        "Subject ID": subject_ids,
        "Percentage of identical matches": [100.0] * len(subject_ids),
        "Alignment length": [139] * len(subject_ids),
        "Number of mismatches": [0] * len(subject_ids),
        "Number of gap openings": [0] * len(subject_ids),
        "Start of alignment in query": [1] * len(subject_ids),
        "End of alignment in query": [139] * len(subject_ids),
        "Start of alignment in subject": [1] * len(subject_ids),
        "End of alignment in subject": [139] * len(subject_ids),
        "Expected value": [1e-100] * len(subject_ids),
        "Bit score": [300] * len(subject_ids),
    })
    return Table(df)


class TestECnumber(BaseTestCase):

    def _run(self, subject_ids: list, taxonomy: str = "fungi") -> dict:
        input_table = _make_diamond_table(subject_ids)
        runner = TaskRunner(
            task_type=ECnumber,
            inputs={"input_table": input_table},
            params={"taxonomy": taxonomy},
        )
        return runner.run()

    # ------------------------------------------------------------------
    # Known fungi proteins → EC numbers should be found
    # ------------------------------------------------------------------

    def test_known_fungi_ids(self):
        """Known S. cerevisiae UniProt IDs → EC numbers found in fungi DB."""
        self.print(f"Test ECnumber: {len(SUBJECT_IDS)} input IDs  ({SUBJECT_IDS})")

        outputs = self._run(SUBJECT_IDS)

        ec_table: Table = outputs["output_table"]
        self.assertIsNotNone(ec_table)
        ec_list = list(ec_table.get_data()["EC number"])
        self.print(f"EC numbers found: {len(ec_list)}")
        self.print("\n" + ec_table.get_data().to_string())
        # At least one real EC number should be found from the 3 known IDs
        self.assertGreater(len(ec_list), 0, msg="Expected at least one EC number hit")

    # ------------------------------------------------------------------
    # Unknown IDs → empty result table
    # ------------------------------------------------------------------

    def test_unknown_ids_returns_empty(self):
        """Completely unknown IDs → output table has 0 rows."""
        fake_ids = ["FAKE001.1", "FAKE002.2", "FAKE003.3"]
        self.print(f"Test ECnumber: unknown IDs {fake_ids}")

        outputs = self._run(fake_ids)

        ec_table: Table = outputs["output_table"]
        self.assertIsNotNone(ec_table)
        rows = ec_table.get_data().shape[0]
        self.print(f"EC numbers found: {rows}  (expected: 0)")
        self.assertEqual(rows, 0, msg="Expected 0 EC hits for unknown IDs")

    # ------------------------------------------------------------------
    # Single known ID → exactly its EC number(s)
    # ------------------------------------------------------------------

    def test_single_known_id(self):
        """Single known ADH1_YEAST ID → EC 1.1.1.1 found."""
        self.print("Test ECnumber: single ID P00330.5 (ADH1_YEAST)")

        outputs = self._run(["P00330.5"])

        ec_table: Table = outputs["output_table"]
        self.assertIsNotNone(ec_table)
        ec_list = list(ec_table.get_data()["EC number"])
        self.print(f"EC numbers: {ec_list}")
        self.assertGreater(len(ec_list), 0, msg="Expected EC number for P00330 (ADH1)")

    # ------------------------------------------------------------------
    # Real Diamond output CSV file as File input
    # ------------------------------------------------------------------

    def test_from_real_diamond_csv_file(self):
        """Real Diamond output CSV (input_ec_diamond_sce.csv) passed as File input."""
        testdata_dir = Settings.get_instance().get_variable("gws_omix", "testdata_dir")
        csv_path = os.path.join(testdata_dir, "alignment", "input_ec_diamond_sce.csv")
        self.print(f"Test ECnumber from File: {csv_path}")

        runner = TaskRunner(
            task_type=ECnumber,
            inputs={"input_table": File(csv_path)},
            params={"taxonomy": "fungi"},
        )
        outputs = runner.run()

        ec_table: Table = outputs["output_table"]
        self.assertIsNotNone(ec_table)
        ec_list = list(ec_table.get_data()["EC number"])
        self.print(f"EC numbers found: {len(ec_list)}")
        self.print("\n" + ec_table.get_data().to_string())
        self.assertGreater(len(ec_list), 0, msg="Expected EC numbers from real Diamond output")

        # Check the text file output — header + 4 EC numbers = 5 lines for this specific input
        ec_file: File = outputs["output_file"]
        self.assertIsNotNone(ec_file)
        with open(ec_file.path) as f:
            lines = [l.strip() for l in f if l.strip()]
        self.print(f"Lines in ec_numbers.txt: {lines}")
        self.assertEqual(len(lines), 5, msg="Expected header + 4 EC numbers = 5 lines for input_ec_diamond_sce.csv")
