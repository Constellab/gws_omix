import os
import re
import tempfile

from gws_core import BaseTestCase, File, TaskRunner,Settings
from gws_omix.kegg.kegg_visualisation import KEGGVisualisation
settings = Settings.get_instance()

class TestKEGGVisualisation(BaseTestCase):

    def _make_no_fc_csv_force_comma(self, src_path: str, header: str) -> File:
        """
        Format:
            NCBI GeneID,_
            1234,
            5678,
        """
        with open(src_path, "r", encoding="utf-8") as f:
            raw_lines = [ln.strip() for ln in f if ln.strip()]

        ids = []
        for ln in raw_lines:
            # skip header-like line if present
            if ln.strip() == header:
                continue
            # accept pure digits only (Entrez IDs)
            if re.fullmatch(r"\d+", ln):
                ids.append(ln)

        fd, out_path = tempfile.mkstemp(suffix=".csv")
        os.close(fd)

        with open(out_path, "w", encoding="utf-8") as out:
            out.write(f"{header},_\n")      # second column exists but is irrelevant
            for gid in ids:
                out.write(f"{gid},\n")      # empty second field

        return File(out_path)

    def _run_kegg(self, infile: File, organism: str, col_entrez: str, foldchange_cols: str):
        runner = TaskRunner(
            task_type=KEGGVisualisation,
            inputs={"deg_file": infile},
            params={
                "organism_name": organism,
                "col_entrez": col_entrez,
                "foldchange_cols": foldchange_cols,
            },
        )
        return runner.run()

    def _print_png_count(self, outputs, label: str):
        pathways = outputs["pathways"].get_resources()
        self.print(f"---- Test: [{label}] PNG generated: {len(pathways)} ----")
        return pathways

    def test_hsa_no_fc(self):
        self.print("Test KEGG Visualisation: HSA, 0FC")

        data_dir = settings.get_variable("gws_omix","testdata_dir")
        src_path = os.path.join(data_dir, "kegg/genes_human.txt")

        # IMPORTANT: make a CSV that forces comma detection
        deg_file = self._make_no_fc_csv_force_comma(src_path, header="NCBI GeneID")

        outputs = self._run_kegg(
            deg_file,
            organism="Homo sapiens",
            col_entrez="NCBI GeneID",
            foldchange_cols="",
        )

        pathways = self._print_png_count(outputs, "HSA/NO_FC")
        self.assertTrue("hsa00480.pathview.png" in pathways)

    def test_hsa_one_fc(self):
        self.print("Test KEGG Visualisation: HSA, 1 FC")
        data_dir = settings.get_variable("gws_omix","testdata_dir")
        deg_file = File(os.path.join(data_dir, "kegg/genes_human_one_fold_change.txt"))

        outputs = self._run_kegg(
            deg_file,
            organism="Homo sapiens",
            col_entrez="NCBI GeneID",
            foldchange_cols="FoldChange",
        )

        pathways = self._print_png_count(outputs, "HSA/1FC")
        self.assertTrue("hsa00480.pathview.png" in pathways)

    def test_hsa_two_fc(self):
        self.print("Test KEGG Visualisation: HSA, 2 FC")
        data_dir = settings.get_variable("gws_omix","testdata_dir")
        deg_file = File(os.path.join(data_dir, "kegg/genes_human_two_fold_change.txt"))

        outputs = self._run_kegg(
            deg_file,
            organism="Homo sapiens",
            col_entrez="NCBI GeneID",
            foldchange_cols="FoldChange,FoldChange 2",
        )

        pathways = self._print_png_count(outputs, "HSA/2FC")
        self.assertTrue("hsa00480.pathview.multi.png" in pathways)