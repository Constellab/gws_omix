import os
import re
import shlex
import urllib.request
from typing import Dict, List, Set, Optional, Tuple

import pandas as pd

from gws_core import (
    ConfigParams,
    ConfigSpecs,
    File,
    InputSpec,
    InputSpecs,
    OutputSpec,
    OutputSpecs,
    ResourceSet,
    ShellProxy,
    StrParam,
    Table,
    TableImporter,
    Task,
    TaskInputs,
    TaskOutputs,
    TypingStyle,
    task_decorator,
)

from .kegg_r_env_task import KeggREnvHelper

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))

# organism autocomplete list
KORG_PATH = os.path.join(_THIS_DIR, "list_organisms_pathview.txt")
_KORG = pd.read_csv(KORG_PATH, sep=None, engine="python", dtype=str, keep_default_na=False)

if "scientific.name" not in _KORG.columns or "kegg.code" not in _KORG.columns:
    raise Exception("list_organisms_pathview.txt must contain columns: scientific.name, kegg.code")

_KORG["scientific.name"] = _KORG["scientific.name"].astype(str).str.strip()
_KORG["kegg.code"] = _KORG["kegg.code"].astype(str).str.strip()

SUPPORTED_ORGANISM_NAMES = sorted([x for x in _KORG["scientific.name"].unique().tolist() if x and x.lower() != "nan"])
NAME_TO_KEGG = dict(zip(_KORG["scientific.name"], _KORG["kegg.code"]))
LOWER_NAME_TO_KEGG = {k.lower(): v for k, v in NAME_TO_KEGG.items() if k and k.lower() != "nan"}


def bh_adjust(pvals: List[float]) -> List[float]:
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    out = [1.0] * m
    prev = 1.0
    for rank in range(m):
        i = order[m - 1 - rank]
        q = pvals[i] * m / (m - rank)
        prev = min(prev, q)
        out[i] = min(1.0, prev)
    return out


def hypergeom_right_tail(N: int, K: int, n: int, k: int) -> float:
    import math
    if k <= 0:
        return 1.0
    max_x = min(K, n)

    def logC(a: int, b: int) -> float:
        return math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)

    denom = logC(N, n)
    p = 0.0
    for x in range(k, max_x + 1):
        num = logC(K, x) + logC(N - K, n - x) - denom
        p += math.exp(num)
    return min(1.0, p)


@task_decorator(
    "KEGGVisualisation",
    human_name="KEGG enrichment + Pathview",
    short_description="Entrez (NCBI) -> KEGG genes -> enrichment -> Pathview (multi-state if multiple FC columns).",
    style=TypingStyle.material_icon(material_icon_name="collections_bookmark", background_color="#d9d9d9"),
)
class KEGGVisualisation(Task):
    """
    Generates a KEGG pathway using the genes specified in the input.

    Please provide a list of genes such as:
    NCBI gene ID,FoldChange
    Gene1,value
    Gene2,value

    If you don't have gene expression, just provide the gene names with a header.
    If you want to compare more than one fold change, add the other fold change in the following columns.
    In the output, each box of the pathway will be separated depending on the number of fold changes.
    So the colours may be different depending on the condition.

    You can find the list of allowed organisms values attached to this story: https://constellab.community/bricks/gws_omix/latest/doc/use-cases/kegg-enrichment-analysis/e601862f-3ad6-47a6-8613-ca378915ca05

    Be aware that this task can take some time, especially the first time, as a virtual environment has to be installed, and also depending on the length of the genes provided, it can take more time.

    In the output you will also get a Table with the pathways where genes are mapped, but these pathways can't be shown with the pathview package.

    KEGG is a database resource for understanding high-level functions and utilities of the biological system.
    Kanehisa Laboratories owns and controls the rights to KEGG.
    Although the KEGG database is made freely available for academic use via the website, it is not in the public domain.
    All commercial use of KEGG requires a license. Please ensure that you have licence to use KEGG database.
    """
    input_specs = InputSpecs({
        "deg_file": InputSpec([File, Table], human_name="DEG table (CSV/TSV)"),
    })

    output_specs = OutputSpecs({
        "pathways": OutputSpec(ResourceSet, human_name="Pathview PNGs"),
        "kegg_enrichment": OutputSpec(Table, human_name="KEGG enrichment (with pathway names)")
    })

    config_specs = ConfigSpecs({
        "organism_name": StrParam(allowed_values=SUPPORTED_ORGANISM_NAMES, human_name="Organism (autocomplete)"),
        "col_entrez": StrParam(
            default_value="NCBI GeneID",
            human_name="NCBI Entrez Gene column",
            short_description="Must contain only digits (NCBI Entrez Gene IDs).",
        ),
        "foldchange_cols": StrParam(
            default_value="",
            human_name="Fold-change columns (optional, ',' separated)",
            short_description="Example: 'log2FoldChange' or 'FoldChange,FoldChange 2'.",
        ),
    })

    # -----------------------------
    def _read_header_only_columns(self, path: str, sep: str) -> List[str]:
        """
        Read ONLY the header (no data rows) to get column names.
        This does not depend on file length.
        """
        df0 = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, nrows=0)
        return [str(c).strip() for c in df0.columns]

    def _detect_sep_and_validate_column_exists(self, path: str, col_entrez: str) -> str:
        """
        Only allowed separators: comma or tab.
        Returns the separator that contains col_entrez in header.
        """
        if not col_entrez:
            raise Exception("Entrez column name is empty.")

        # Try comma, then tab. No sniffing.
        for sep in [",", "\t"]:
            try:
                cols = self._read_header_only_columns(path, sep=sep)
                if col_entrez in cols:
                    return sep
            except Exception:
                continue

        # Build a debug message with what we *did* parse (best-effort)
        parsed_cols: Optional[List[str]] = None
        parsed_sep: Optional[str] = None
        for sep in [",", "\t"]:
            try:
                parsed_cols = self._read_header_only_columns(path, sep=sep)
                parsed_sep = sep
                break
            except Exception:
                pass

        if parsed_cols is not None:
            raise Exception(
                f"Entrez column '{col_entrez}' not found in input. "
                f"Allowed separators are ',' or tab. "
                f"Parsed columns (sep={'TAB' if parsed_sep == chr(9) else ','}): {parsed_cols}"
            )
        raise Exception(
            f"Entrez column '{col_entrez}' not found in input. "
            f"Allowed separators are ',' or tab. Also failed to parse header."
        )

    def _stream_validate_entrez_digits(
        self,
        path: str,
        sep: str,
        col_entrez: str,
        chunksize: int = 1_000_000,
    ) -> None:
        """
        Streaming validation: reads only the Entrez column in chunks.
        Supports extremely large files without RAM blowup.
        Raises early if too many non-numeric values are detected.
        """
        digit_re = re.compile(r"^\d+(\.0)?$")

        bad_count = 0
        nonempty_seen = 0

        # read only the one column
        for chunk in pd.read_csv(
            path,
            sep=sep,
            dtype=str,
            keep_default_na=False,
            usecols=[col_entrez],
            chunksize=chunksize,
        ):
            s = chunk[col_entrez].astype(str).str.strip()
            s = s[s != ""]
            if s.empty:
                continue

            nonempty_seen += int(s.shape[0])

            # vectorized-ish: still uses python for regex per element, but chunked
            bad_mask = ~s.map(lambda v: bool(digit_re.fullmatch(v)))
            bad_here = int(bad_mask.sum())
            bad_count += bad_here

            # Keep your original spirit: if there are many bads, fail early
            if bad_count > 5 and nonempty_seen > 0:
                # show a few examples
                examples = s[bad_mask].head(10).tolist()
                raise Exception(
                    "The selected Entrez column contains non-numeric values.\n"
                    "It must contain NCBI Entrez Gene IDs (digits only).\n"
                    f"Examples: {examples}\n"
                    "Please run OmiX 'Gene ID conversion' with target_namespace='ENTREZGENE' "
                    "(and numeric_namespace='ENTREZGENE_ACC' if your IDs are numeric), "
                    "then use the produced Entrez column."
                )

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        shell: ShellProxy = KeggREnvHelper.create_proxy(self.message_dispatcher)

        organism_name = str(params["organism_name"]).strip()
        specie = self._specie_from_name(organism_name)
        self.log_info_message(f"Selected organism: {organism_name} -> KEGG code: {specie}")

        col_entrez = (params.get("col_entrez") or "").strip()
        foldchange_cols = (params.get("foldchange_cols") or "").strip()

        deg_path = self._materialize_input(inputs["deg_file"], shell)

        # -----------------------------
        # FIXED VALIDATION (no nrows limit, no sep sniffing)
        # -----------------------------
        sep_used = self._detect_sep_and_validate_column_exists(deg_path, col_entrez)
        self._stream_validate_entrez_digits(deg_path, sep_used, col_entrez)

        # 1) Prepare KEGG gene table for pathview (and mapping table)
        prep_py = os.path.join(_THIS_DIR, "_kegg_prepare_kegg_gene_data.py")
        out_gene_kegg = os.path.join(shell.working_dir, "gene_kegg.csv")
        out_used = os.path.join(shell.working_dir, "gene_used.csv")

        cmd = " ".join([
            "python3", shlex.quote(prep_py),
            "--deg", shlex.quote(deg_path),
            "--specie", shlex.quote(specie),
            "--col_entrez", shlex.quote(col_entrez),
            "--foldchange_cols", shlex.quote(foldchange_cols),
            "--out_gene_kegg", shlex.quote(out_gene_kegg),
            "--out_used", shlex.quote(out_used),
        ])

        rc = shell.run(cmd, shell_mode=True)
        if rc != 0 or (not os.path.exists(out_gene_kegg)) or os.path.getsize(out_gene_kegg) == 0:
            raise Exception("Gene preparation failed: gene_kegg.csv not produced.")

        used_df = pd.read_csv(out_used, dtype=str, keep_default_na=False)
        mapped_full = used_df["kegg_full"].astype(str)
        query_full = sorted(set([x for x in mapped_full.tolist() if x and x.startswith(specie + ":")]))

        # 2) Universe + pathway gene sets
        gp = self._kegg_link_pathway_species(specie)
        if gp.empty:
            raise Exception("KEGG link/pathway/<specie> returned empty universe (network/proxy or parsing issue).")

        universe = set(gp["gene"].tolist())
        query = [g for g in query_full if g in universe]
        self.log_info_message(f"Universe KEGG gene size: {len(universe)} | Query∩Universe: {len(query)}")
        if len(query) == 0:
            raise Exception("Query genes do not intersect KEGG universe (after Entrez->KEGG conversion).")

        pw2genes: Dict[str, Set[str]] = {}
        for _, r in gp.iterrows():
            pw2genes.setdefault(r["pathway"], set()).add(r["gene"])

        # 3) Pathway names
        pw_names_df = self._kegg_list_pathways(specie)
        pw_name_map = dict(zip(pw_names_df["pathway"], pw_names_df["pathway_name"]))

        # 4) Enrichment
        qset = set(query)
        N = len(universe)
        n = len(query)

        rows = []
        for pw, geneset in pw2genes.items():
            K = len(geneset)
            k = len(qset & geneset)  # overlap
            if k == 0:
                continue
            p = hypergeom_right_tail(N, K, n, k)
            rows.append((pw, k, K, n, N, p))

        if not rows:
            raise Exception("Enrichment table is empty (no overlaps).")

        enr = pd.DataFrame(rows, columns=[
            "pathway", "overlap", "pathway_size", "query_size", "universe_size", "pvalue"
        ])
        enr["padj"] = bh_adjust(enr["pvalue"].tolist())
        enr["pathway_name"] = enr["pathway"].map(lambda p: pw_name_map.get(p, ""))
        enr = enr.sort_values(["padj", "pvalue", "overlap"], ascending=[True, True, False]).reset_index(drop=True)

        enr_csv = os.path.join(shell.working_dir, "kegg_enrichment.csv")
        enr.to_csv(enr_csv, index=False)

        # 5) Pathview: render ALL overlapped pathways (no top-N limit)
        def to_pid5(pw: str) -> str:
            digs = "".join(ch for ch in str(pw) if ch.isdigit())
            return digs[-5:] if len(digs) >= 5 else ""

        pids = [to_pid5(x) for x in enr["pathway"].astype(str).tolist()]
        pids = [x for x in pids if x]

        pathway_list = os.path.join(shell.working_dir, "pathway_kegg.txt")
        with open(pathway_list, "w") as f:
            f.write("\n".join(pids))

        kegg_dir = os.path.join(shell.working_dir, "kegg_cache")
        os.makedirs(kegg_dir, exist_ok=True)
        self._prefetch_kegg_files(specie, pids, kegg_dir)

        r_script = os.path.join(_THIS_DIR, "kegg_visualisation.R")

        cmd_r = " ".join([
            "bash", "-c",
            shlex.quote(
                f"cd {shlex.quote(shell.working_dir)} && "
                f"Rscript --vanilla {shlex.quote(r_script)} "
                f"{shlex.quote(out_gene_kegg)} {shlex.quote(specie)} "
                f"{shlex.quote(pathway_list)} {shlex.quote(kegg_dir)}"
            )
        ])

        rc = shell.run(cmd_r, shell_mode=True)
        if rc != 0:
            raise Exception("R/pathview crashed (non-zero exit). Check list_pathway_error.csv.")

        # collect pngs
        rs = ResourceSet()
        rs.name = "Pathview images"
        for fn in os.listdir(shell.working_dir):
            fp = os.path.join(shell.working_dir, fn)
            if os.path.isfile(fp) and (fn.endswith(".pathview.png") or fn.endswith(".pathview.multi.png")):
                rs.add_resource(File(fp), fn)

        # import tables
        kegg_enrichment_tbl = TableImporter.call(File(enr_csv), params={"index_column": -1})

        err_fp = os.path.join(shell.working_dir, "list_pathway_error.csv")
        if os.path.exists(err_fp) and os.path.getsize(err_fp) > 0:
            _ = TableImporter.call(File(err_fp), params={"index_column": -1})
        else:
            pd.DataFrame(columns=["pathway", "status", "message"]).to_csv(err_fp, index=False)

        return {
            "pathways": rs,
            "kegg_enrichment": kegg_enrichment_tbl,
        }

    def _specie_from_name(self, organism_name: str) -> str:
        if organism_name in NAME_TO_KEGG:
            return NAME_TO_KEGG[organism_name]
        low = organism_name.lower()
        if low in LOWER_NAME_TO_KEGG:
            return LOWER_NAME_TO_KEGG[low]
        raise Exception(f"Organism '{organism_name}' not found in list_organisms_pathview.txt")

    def _materialize_input(self, inp, shell: ShellProxy) -> str:
        if isinstance(inp, Table):
            p = os.path.join(shell.working_dir, "deg_input.csv")
            inp.to_dataframe().to_csv(p, index=False)
            return p
        return inp.path

    def _fetch(self, url: str, timeout: int = 180) -> str:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            return resp.read().decode("utf-8", errors="replace")

    def _kegg_link_pathway_species(self, specie: str) -> pd.DataFrame:
        url = f"https://rest.kegg.jp/link/pathway/{specie}"
        raw = self._fetch(url, timeout=240)

        head = raw.lstrip()[:300].lower()
        if head.startswith("<!doctype") or head.startswith("<html") or "forbidden" in head or "error" in head:
            sample = "\n".join(raw.splitlines()[:10])
            raise Exception(f"Unexpected KEGG response for {url}. First lines:\n{sample}")

        genes, pws = [], []
        for line in raw.splitlines():
            line = line.strip()
            if not line:
                continue

            if "\t" in line:
                left, right = line.split("\t", 1)
            else:
                parts = line.split(None, 1)
                if len(parts) < 2:
                    continue
                left, right = parts[0], parts[1]

            left, right = left.strip(), right.strip()

            if left.startswith(specie + ":") and right.startswith("path:"):
                gene = left
                pw = right.replace("path:", "")
            elif right.startswith(specie + ":") and left.startswith("path:"):
                gene = right
                pw = left.replace("path:", "")
            else:
                continue

            if re.fullmatch(rf"{re.escape(specie)}\d{{5}}", pw):
                genes.append(gene)
                pws.append(pw)

        return pd.DataFrame({"gene": genes, "pathway": pws})

    def _kegg_list_pathways(self, specie: str) -> pd.DataFrame:
        url = f"https://rest.kegg.jp/list/pathway/{specie}"
        raw = self._fetch(url, timeout=240)

        ids, names = [], []
        for line in raw.splitlines():
            line = line.strip()
            if not line:
                continue
            a, b = line.split("\t", 1)
            pid = a.replace("path:", "").strip()
            name = b.strip()
            if re.fullmatch(rf"{re.escape(specie)}\d{{5}}", pid):
                ids.append(pid)
                names.append(name)

        return pd.DataFrame({"pathway": ids, "pathway_name": names})

    def _download(self, url: str, out_path: str) -> None:
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:
                data = resp.read()
            with open(out_path, "wb") as f:
                f.write(data)
        except Exception:
            pass

    def _prefetch_kegg_files(self, specie: str, pathways_5: List[str], kegg_dir: str) -> None:
        for pid in pathways_5:
            base = f"{specie}{pid}"
            xml_path = os.path.join(kegg_dir, f"{base}.xml")
            png_path = os.path.join(kegg_dir, f"{base}.png")
            if not os.path.exists(xml_path):
                self._download(f"https://rest.kegg.jp/get/{base}/kgml", xml_path)
            if not os.path.exists(png_path):
                self._download(f"https://rest.kegg.jp/get/{base}/image", png_path)