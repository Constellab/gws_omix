#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
_gsea.py — Pre-ranked GSEA (GMT) using gseapy.

Inputs:
  - DE results table (CSV/TSV-like; sep auto-detected)
  - genes_colname: gene identifiers
  - stat_colname: signed ranking statistic (e.g., Wald stat / t-stat)
  - gmt_file: ONE GMT file

Outputs (written into csv_dir / fig_dir):
  - csv/GSEA_results.csv
  - csv/GSEA_matrix_long.csv
  - figs/GSEA_top20_barplot.png
  - figs/GSEA_enrichment_01_<TERM>.png ... figs/GSEA_enrichment_XX_<TERM>.png

NEW:
  - --n_enrichplots : number of classic enrichment plots to generate (default=10)
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import gseapy as gp


# ── fixed knobs (not exposed) ─────────────────────────────────────
SEED: int = 0
SCORE_TYPE: str = "std"
GSEA_WEIGHT: float = 1.0
TOPN_BARPLOT: int = 20
DEFAULT_N_CLASSIC_PLOTS: int = 10  # <- exposed via CLI

SAFE_NAME_RE = re.compile(r"[^A-Za-z0-9_.-]+")


def eprint(*args) -> None:
    print(*args, file=sys.stderr)


def safe_name(s: str, max_len: int = 90) -> str:
    s = (s or "").strip().replace(" ", "_")
    s = SAFE_NAME_RE.sub("_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s[:max_len] if len(s) > max_len else s


def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg adjusted p-values (padj)."""
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    finite = np.isfinite(p)
    if not finite.any():
        return out

    p2 = p[finite]
    m = len(p2)
    order = np.argsort(p2, kind="mergesort")
    ranks = np.arange(1, m + 1, dtype=float)

    q = p2[order] * (m / ranks)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)

    out_idx = np.where(finite)[0][order]
    out[out_idx] = q
    return out


def collapse_duplicates_max_abs(df: pd.DataFrame, gene_col: str, stat_col: str) -> pd.DataFrame:
    d = df[[gene_col, stat_col]].copy()
    d[gene_col] = d[gene_col].astype(str)
    d[stat_col] = pd.to_numeric(d[stat_col], errors="coerce")
    d = d.dropna(subset=[gene_col, stat_col])
    if d.empty:
        return pd.DataFrame(columns=["gene", "stat"])

    d["_abs_"] = d[stat_col].abs()
    idx = d.groupby(gene_col)["_abs_"].idxmax()
    out = d.loc[idx, [gene_col, stat_col]].copy()
    out = out.rename(columns={gene_col: "gene", stat_col: "stat"})
    out = out.sort_values("stat", ascending=False).reset_index(drop=True)
    return out


def read_gmt_meta(gmt_path: Path) -> Tuple[Dict[str, str], Dict[str, int], set[str], int]:
    term_desc: Dict[str, str] = {}
    term_size: Dict[str, int] = {}
    universe: set[str] = set()
    n_sets = 0

    with gmt_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            term = parts[0].strip()
            desc = parts[1].strip() if len(parts) > 1 else ""
            genes = [g.strip() for g in parts[2:] if g.strip()]
            if not term or not genes:
                continue
            n_sets += 1
            term_desc[term] = desc
            term_size[term] = len(genes)
            universe.update(genes)

    return term_desc, term_size, universe, n_sets


def normalize_gseapy_table(res2d: pd.DataFrame) -> pd.DataFrame:
    d = res2d.copy()

    if "Term" not in d.columns:
        d = d.reset_index().rename(columns={"index": "set_id"})
    else:
        d = d.rename(columns={"Term": "set_id"})

    rename = {}
    for c in list(d.columns):
        lc = str(c).lower().strip()
        if lc == "es":
            rename[c] = "ES"
        elif lc == "nes":
            rename[c] = "NES"
        elif lc in {"pval", "pvalue", "nom p-val", "nom_pval", "nominal_p_value"}:
            rename[c] = "p_value"
        elif lc in {"fdr", "fdr q-val", "fdr_qval", "qval", "q_value", "fdr q val"}:
            rename[c] = "qvalue"
        elif lc in {"lead_genes", "leading_edge"}:
            rename[c] = "lead_genes"
        elif lc in {"geneset_size", "size"}:
            rename[c] = "set_size"

    d = d.rename(columns=rename)

    for col in ["set_id", "ES", "NES", "p_value", "qvalue", "lead_genes"]:
        if col not in d.columns:
            d[col] = np.nan

    d["ES"] = pd.to_numeric(d["ES"], errors="coerce")
    d["NES"] = pd.to_numeric(d["NES"], errors="coerce")
    d["p_value"] = pd.to_numeric(d["p_value"], errors="coerce")
    d["qvalue"] = pd.to_numeric(d["qvalue"], errors="coerce")
    if "set_size" in d.columns:
        d["set_size"] = pd.to_numeric(d["set_size"], errors="coerce")

    rankp = d["qvalue"].copy()
    if rankp.isna().all():
        rankp = d["p_value"].copy()

    d["_rank_p_"] = rankp
    d["_abs_nes"] = d["NES"].abs()
    d = d.dropna(subset=["_rank_p_"]).sort_values(["_rank_p_", "_abs_nes"], ascending=[True, False])
    d = d.drop(columns=["_rank_p_", "_abs_nes"])
    return d


def leading_edge_size_and_ratio(lead_genes_val, set_size: int) -> Tuple[int, float]:
    if set_size <= 0:
        return 0, 0.0
    if lead_genes_val is None or (isinstance(lead_genes_val, float) and np.isnan(lead_genes_val)):
        return 0, 0.0

    s = str(lead_genes_val).strip()
    if not s:
        return 0, 0.0

    if "/" in s:
        genes = [g.strip() for g in s.split("/") if g.strip()]
    else:
        genes = [g.strip() for g in re.split(r"[\s,;]+", s) if g.strip()]

    n = len(genes)
    return n, float(n / set_size)


def make_barplot(df: pd.DataFrame, outpng: Path, topn: int = TOPN_BARPLOT) -> None:
    d = df.copy()
    d = d.dropna(subset=["NES"])
    if d.empty:
        return

    rankp = pd.to_numeric(d.get("qvalue", np.nan), errors="coerce")
    if np.all(~np.isfinite(rankp)):
        rankp = pd.to_numeric(d.get("padj", np.nan), errors="coerce")
    if np.all(~np.isfinite(rankp)):
        rankp = pd.to_numeric(d.get("p_value", np.nan), errors="coerce")

    d["_rank_p_"] = rankp
    d["_abs_nes"] = d["NES"].abs()
    d = d.dropna(subset=["_rank_p_"]).sort_values(["_rank_p_", "_abs_nes"], ascending=[True, False]).head(topn)
    d = d.iloc[::-1]

    plt.figure(figsize=(10.0, max(3.2, 0.40 * len(d))))
    y = np.arange(len(d))
    plt.barh(y, d["NES"].values)
    plt.yticks(y, d["set_id"].astype(str).values)
    plt.xlabel("NES")
    plt.title(f"GSEA top {min(topn, len(d))} gene sets")
    plt.tight_layout()
    plt.savefig(outpng, dpi=300)
    plt.close()


def build_long_matrix(res_df: pd.DataFrame, ranked: pd.DataFrame) -> pd.DataFrame:
    cols = ["set_id", "gene_label", "NES", "p_value", "padj", "qvalue", "stat"]
    if res_df is None or res_df.empty:
        return pd.DataFrame(columns=cols)

    stat_map = dict(zip(ranked["gene"].astype(str), ranked["stat"].astype(float)))

    def parse_le(x) -> List[str]:
        if x is None or (isinstance(x, float) and np.isnan(x)):
            return []
        s = str(x).strip()
        if not s:
            return []
        if "/" in s:
            return [g.strip() for g in s.split("/") if g.strip()]
        return [g.strip() for g in re.split(r"[\s,;]+", s) if g.strip()]

    rows = []
    for _, r in res_df.iterrows():
        le = parse_le(r.get("lead_genes", None))
        for g in le:
            rows.append({
                "set_id": str(r["set_id"]),
                "gene_label": g,
                "NES": float(r.get("NES", np.nan)),
                "p_value": float(r.get("p_value", np.nan)),
                "padj": float(r.get("padj", np.nan)),
                "qvalue": float(r.get("qvalue", np.nan)),
                "stat": float(stat_map.get(g, np.nan)),
            })
    return pd.DataFrame(rows, columns=cols)


def _plot_gsea_classic(pre: Any, term_key: Any, outpng: Path) -> None:
    """
    gseapy plot API differs by version. We support both signatures.
    """
    # Try "new" signature first
    try:
        gp.plot.gseaplot2(
            term=term_key,
            rank_metric=pre.ranking,
            **pre.results[term_key],
            ofname=str(outpng),
        )
        return
    except TypeError:
        pass
    except Exception:
        pass

    # Old signature: terms + RESs
    try:
        gp.plot.gseaplot2(
            terms=[term_key],
            RESs=[pre.results[term_key]],
            rank_metric=pre.ranking,
            ofname=str(outpng),
        )
        return
    except TypeError:
        pass

    # Fallback
    gp.plot.gseaplot(
        term=term_key,
        rank_metric=pre.ranking,
        **pre.results[term_key],
        ofname=str(outpng),
    )


def main() -> int:
    ap = argparse.ArgumentParser(description="Pre-ranked GSEA (single GMT) using gseapy.")
    ap.add_argument("--infile", required=True)
    ap.add_argument("--genes_colname", required=True)
    ap.add_argument("--stat_colname", required=True)
    ap.add_argument("--gmt_file", required=True)

    ap.add_argument("--outdir", default=".")
    ap.add_argument("--csv_dir", default=None)
    ap.add_argument("--fig_dir", default=None)

    ap.add_argument("--min_size", default="15")
    ap.add_argument("--max_size", default="500")
    ap.add_argument("--nperm", default="10000")
    ap.add_argument("--processes", default="4")

    # ✅ NEW
    ap.add_argument("--n_enrichplots", default=str(DEFAULT_N_CLASSIC_PLOTS))

    args = ap.parse_args()

    outdir = Path(args.outdir)
    csv_dir = Path(args.csv_dir) if args.csv_dir else (outdir / "csv")
    fig_dir = Path(args.fig_dir) if args.fig_dir else (outdir / "figs")
    csv_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    infile = Path(args.infile)
    gmt_file = Path(args.gmt_file)
    if not infile.exists():
        eprint(f"ERROR: infile not found: {infile}")
        return 2
    if not gmt_file.exists():
        eprint(f"ERROR: gmt_file not found: {gmt_file}")
        return 2

    # auto-sep read
    df = pd.read_csv(str(infile), sep=None, engine="python")
    gene_col = args.genes_colname
    stat_col = args.stat_colname
    if gene_col not in df.columns:
        eprint(f"ERROR: genes_colname '{gene_col}' not found in input.")
        return 2
    if stat_col not in df.columns:
        eprint(f"ERROR: stat_colname '{stat_col}' not found in input.")
        return 2

    ranked = collapse_duplicates_max_abs(df, gene_col, stat_col)
    if ranked.empty:
        eprint("ERROR: No valid (gene, stat) after cleaning.")
        return 3

    term_desc, term_size, gmt_universe, n_sets = read_gmt_meta(gmt_file)
    ranked_genes = set(ranked["gene"].astype(str).tolist())
    overlap = len(ranked_genes.intersection(gmt_universe))
    eprint(f"[INFO] ranked_genes={len(ranked_genes)} gmt_sets={n_sets} gmt_unique_genes={len(gmt_universe)} overlap_genes={overlap}")
    if overlap == 0:
        eprint("[WARN] Zero overlap between ranked genes and GMT gene universe. Check ID namespace.")

    gseapy_work = outdir / "_gseapy_run"
    gseapy_work.mkdir(parents=True, exist_ok=True)

    min_size = int(float(args.min_size))
    max_size = int(float(args.max_size))
    nperm = int(float(args.nperm))
    processes = int(float(args.processes))

    try:
        n_classic_plots = int(float(args.n_enrichplots))
    except Exception:
        n_classic_plots = DEFAULT_N_CLASSIC_PLOTS
    n_classic_plots = max(0, n_classic_plots)

    try:
        pre = gp.prerank(
            rnk=ranked[["gene", "stat"]],
            gene_sets=str(gmt_file),
            processes=max(1, processes),
            permutation_num=nperm,
            outdir=str(gseapy_work),
            seed=SEED,
            min_size=min_size,
            max_size=max_size,
            ascending=False,
            verbose=False,
            weighted_score_type=GSEA_WEIGHT,
            score_type=SCORE_TYPE,
            format="png",
        )
    except Exception as ex:
        eprint(f"ERROR: gseapy.prerank failed: {ex}")
        return 4

    res2d = getattr(pre, "res2d", None)
    if res2d is None or res2d.empty:
        eprint("[WARN] No results produced by gseapy (check overlap + min/max size).")
        (csv_dir / "GSEA_results.csv").write_text("", encoding="utf-8")
        (csv_dir / "GSEA_matrix_long.csv").write_text("", encoding="utf-8")
        return 0

    core = normalize_gseapy_table(res2d)

    # Add GMT description + authoritative set size
    core["description"] = core["set_id"].map(lambda t: term_desc.get(str(t), ""))
    core["set_size"] = core["set_id"].map(lambda t: int(term_size.get(str(t), np.nan)) if str(t) in term_size else np.nan)

    # padj from nominal p_value (always)
    core["padj"] = bh_adjust(pd.to_numeric(core["p_value"], errors="coerce").to_numpy())

    # qvalue from gseapy if present; else fallback to padj
    if pd.to_numeric(core["qvalue"], errors="coerce").isna().all():
        core["qvalue"] = core["padj"]

    # Leading edge metrics
    le_sizes, le_ratios = [], []
    for _, r in core.iterrows():
        ss = r.get("set_size", np.nan)
        ss_int = int(ss) if np.isfinite(ss) else 0
        n, ratio = leading_edge_size_and_ratio(r.get("lead_genes", None), ss_int)
        le_sizes.append(n)
        le_ratios.append(ratio)
    core["leading_edge_size"] = le_sizes
    core["leading_edge_ratio"] = le_ratios

    out_cols = [
        "set_id", "description", "set_size",
        "ES", "NES", "p_value", "padj", "qvalue",
        "lead_genes", "leading_edge_size", "leading_edge_ratio",
    ]
    for c in out_cols:
        if c not in core.columns:
            core[c] = np.nan
    out_table = core[out_cols].copy()
    out_table.to_csv(csv_dir / "GSEA_results.csv", index=False)

    long_df = build_long_matrix(out_table, ranked)
    if long_df.empty:
        (csv_dir / "GSEA_matrix_long.csv").write_text("", encoding="utf-8")
    else:
        long_df.drop_duplicates(subset=["set_id", "gene_label"]).to_csv(csv_dir / "GSEA_matrix_long.csv", index=False)

    # Barplot
    make_barplot(out_table, fig_dir / "GSEA_top20_barplot.png", topn=TOPN_BARPLOT)

    # Classic plots (✅ now uses n_classic_plots)
    results_keys = list(getattr(pre, "results", {}).keys())
    key_map = {str(k).strip(): k for k in results_keys}

    d_plot = out_table.dropna(subset=["NES"]).copy()
    rankp = pd.to_numeric(d_plot["qvalue"], errors="coerce")
    if rankp.isna().all():
        rankp = pd.to_numeric(d_plot["padj"], errors="coerce")
    if rankp.isna().all():
        rankp = pd.to_numeric(d_plot["p_value"], errors="coerce")

    d_plot["_rank_p_"] = rankp
    d_plot["_abs_nes"] = d_plot["NES"].abs()
    d_plot = (
        d_plot.dropna(subset=["_rank_p_"])
        .sort_values(["_rank_p_", "_abs_nes"], ascending=[True, False])
        .head(n_classic_plots)
    )

    for i, term in enumerate(d_plot["set_id"].astype(str).tolist(), start=1):
        lookup = key_map.get(term.strip())
        if lookup is None and term in getattr(pre, "results", {}):
            lookup = term
        if lookup is None:
            eprint(f"[WARN] Cannot find term key for plotting: {term}")
            continue

        term_safe = safe_name(str(term), max_len=70)
        outp = fig_dir / f"GSEA_enrichment_{i:02d}_{term_safe}.png"
        try:
            _plot_gsea_classic(pre, lookup, outp)
        except Exception as ex:
            eprint(f"[WARN] Failed to plot term='{term}': {ex}")

    eprint("✔ Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())