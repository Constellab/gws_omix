from gws_core import BaseTestCase, TaskRunner

from gws_omix.file.fastq_download_sra_ena.download_sra_ena import FastqDLRunner

# ---------------------------------------------------------------------------
# Accession to download
# ---------------------------------------------------------------------------
ACCESSION = "SRR7739756"


class TestFastqDLRunner(BaseTestCase):

    def _run(self, accession: str, provider: str = "sra", cpus: int = 2):
        runner = TaskRunner(
            task_type=FastqDLRunner,
            inputs={},
            params={
                "accession": accession,
                "provider": provider,
                "cpus": cpus,
            },
        )
        return runner.run()

    # ------------------------------------------------------------------
    # SRA provider
    # ------------------------------------------------------------------

    def test_sra(self):
        """Download accession from SRA."""
        self.print(f"Test FastqDLRunner: {ACCESSION} via SRA")

        outputs = self._run(ACCESSION, provider="sra")

        fastq_dir = outputs["fastq_dir"]
        run_info = outputs["run_info_tbl"]
        self.assertIsNotNone(fastq_dir)
        self.assertIsNotNone(run_info)
        self.assertGreaterEqual(run_info.get_data().shape[0], 1)
        self.print(f"FASTQ dir: {fastq_dir.path}")
        self.print(f"Run info columns: {list(run_info.get_data().columns)}")

    # ------------------------------------------------------------------
    # ENA provider
    # ------------------------------------------------------------------

    def test_ena(self):
        """Download accession from ENA."""
        self.print(f"Test FastqDLRunner: {ACCESSION} via ENA")

        outputs = self._run(ACCESSION, provider="ena")

        self.assertIsNotNone(outputs["fastq_dir"])
        self.assertIsNotNone(outputs["run_info_tbl"])
        self.assertGreaterEqual(outputs["run_info_tbl"].get_data().shape[0], 1)
