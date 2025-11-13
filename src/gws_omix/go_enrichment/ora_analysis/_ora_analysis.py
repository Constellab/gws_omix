#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
import argparse, ast, os, re
from pathlib import Path
from typing import Sequence, Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
import seaborn as sns
from gprofiler import GProfiler

# ─────────────────────────── Helpers ────────────────────────────
ENSEMBL_RE = r"^ENS[A-Z]*G\d+(\.\d+)?$"

def guess_code(scientific_name: str) -> str:
    parts = re.split(r"[\s_]+", scientific_name.strip())
    if len(parts) < 2:
        raise ValueError(f"Scientific name must be 'Genus species', got: '{scientific_name}'")
    return (parts[0][0] + parts[1]).lower()

def strip_version(ids: Sequence[str]) -> list[str]:
    return [re.sub(r"\.\d+$", "", str(x)) for x in ids]

def looks_like_ensembl(s: pd.Series, thr: float = 0.3) -> bool:
    if s is None or s.empty:
        return False
    return s.astype(str).str.match(ENSEMBL_RE).mean() > thr

def uniq(seq: Sequence[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for x in seq:
        if x not in seen:
            seen.add(x); out.append(x)
    return out

def series_clean_ids(s: pd.Series) -> pd.Series:
    s = s.dropna().astype(str)
    if looks_like_ensembl(s, thr=0.3):
        s = pd.Series(strip_version(s), index=s.index)
    return s

def _score_and_pval_cols(df: pd.DataFrame) -> Tuple[pd.Series, pd.Series]:
    """Return score (-log10 p) and p-value series from g:Profiler results."""
    if "adjusted_p_value" in df.columns:
        adjp = pd.to_numeric(df["adjusted_p_value"], errors="coerce")
        score = -np.log10(adjp.clip(lower=1e-300))
        pval  = pd.to_numeric(df.get("p_value", adjp), errors="coerce")
    elif "p_value" in df.columns:
        pval  = pd.to_numeric(df["p_value"], errors="coerce")
        score = -np.log10(pval.clip(lower=1e-300))
    elif "negative_log10_of_adjusted_p_value" in df.columns:
        score = pd.to_numeric(df["negative_log10_of_adjusted_p_value"], errors="coerce")
        pval  = 10 ** (-score)
    else:
        pval = pd.Series(np.nan, index=df.index)
        score = pd.Series(np.nan, index=df.index)
    return score, pval

def _intersections_column(df: pd.DataFrame) -> str | None:
    # more permissive to accommodate lib / API variants
    for c in ("intersections", "intersecting_genes", "queryGenes",
              "genes", "overlap", "overlapping_genes"):
        if c in df.columns:
            return c
    return None

# ───────────────────── Static barplot (optionnel) ─────────────────
def barplot(df: pd.DataFrame, outpng: Path, title: str, topn=20) -> bool:
    if df is None or df.empty:
        return False
    score, _ = _score_and_pval_cols(df)
    if score.isna().all():
        return False
    d = df.copy()
    d["_score_"] = score
    labels = (d["name"] if "name" in d.columns else d.get("term")).astype(str)
    top = d.assign(_label_=labels).sort_values("_score_", ascending=False).head(topn)
    if top.empty:
        return False
    plt.figure(figsize=(7.8, max(2.5, 0.45 * len(top))))
    sns.barplot(data=top, x="_score_", y="_label_", hue="source", dodge=False)
    plt.xlabel("-log10(adj p)")
    plt.ylabel("Term")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpng, dpi=300)
    plt.close()
    return True

# ───────────────────── g:Profiler & filtering ────────────────────
def run_gprofiler(query_genes: list[str], organism_code: str, sources: list[str]) -> pd.DataFrame:
    if not query_genes:
        return pd.DataFrame()
    gp = GProfiler(return_dataframe=True)
    res = gp.profile(
        organism=organism_code,
        query=query_genes,
        all_results=True,
        sources=sources,
        user_threshold=1.0,
        no_evidences=False,
        ordered=False,
        background=None,
    )
    return res if res is not None else pd.DataFrame()

def filter_like_website(res: pd.DataFrame, alpha: float) -> pd.DataFrame:
    if res is None or res.empty:
        return res
    r = res.copy()
    if "adjusted_p_value" in r.columns:
        r = r[r["adjusted_p_value"] <= alpha]
    elif "p_value" in r.columns:
        r = r[r["p_value"] <= alpha]
    elif "negative_log10_of_adjusted_p_value" in r.columns:
        r = r[r["negative_log10_of_adjusted_p_value"] >= -np.log10(alpha)]
    else:
        r = pd.DataFrame(columns=r.columns)
    return r

# ───────────────────── Long matrix (avec log2fc) ─────────────────
def explode_long(df: pd.DataFrame, de_lfc_map: Dict[str, float]) -> pd.DataFrame:
    cols = ["term_name", "term_id", "source", "gene_label", "score", "pvalue", "log2fc"]
    if df is None or df.empty:
        return pd.DataFrame(columns=cols)

    score, pval = _score_and_pval_cols(df)
    inter_col = _intersections_column(df)
    if inter_col is None:
        return pd.DataFrame(columns=cols)

    def _to_list(x):
        if isinstance(x, (list, tuple, set)): return list(x)
        if isinstance(x, str):
            xs = x.strip()
            if xs.startswith("[") or xs.startswith("("):
                try:
                    return list(ast.literal_eval(xs))
                except Exception:
                    pass
            return re.split(r"[\s,;]+", xs.strip())
        return []

    rows = []
    for i, genes in df[inter_col].apply(_to_list).items():
        if not genes:
            continue
        term = df.at[i, "name"] if "name" in df.columns else df.at[i, "term"] if "term" in df.columns else ""
        src  = df.at[i, "source"] if "source" in df.columns else ""
        tid  = df.at[i, "native"] if "native" in df.columns else (df.at[i, "term_id"] if "term_id" in df.columns else "")
        sc   = float(score.loc[i]) if i in score.index else np.nan
        pv   = float(pval.loc[i])  if i in pval.index  else np.nan
        for g in genes:
            gstr = str(g).strip()
            lfc  = float(de_lfc_map.get(gstr, np.nan))
            rows.append({
                "term_name": str(term), "term_id": str(tid), "source": str(src),
                "gene_label": gstr, "score": sc, "pvalue": pv, "log2fc": lfc
            })
    return pd.DataFrame(rows)

# ─────────────── Build ALL / UP / DOWN gene sets ────────────────
def build_gene_sets(df: pd.DataFrame, id_series: pd.Series, padj_thr: float, lfc_thr: float) -> Dict[str, List[str]]:
    ids = series_clean_ids(id_series)
    if "padj" in df.columns:
        padj = pd.to_numeric(df["padj"], errors="coerce")
        mask_all = padj <= padj_thr
    elif "pvalue" in df.columns:
        pval = pd.to_numeric(df["pvalue"], errors="coerce")
        mask_all = pval <= padj_thr
    else:
        mask_all = pd.Series(True, index=df.index)

    result: Dict[str, List[str]] = {}
    result["ALL_DE"] = uniq(ids.loc[mask_all[mask_all].index].tolist())

    if "log2FoldChange" in df.columns:
        lfc = pd.to_numeric(df["log2FoldChange"], errors="coerce")
        up_mask = mask_all & (lfc >= max(0.0, lfc_thr))
        dn_mask = mask_all & (lfc <= -max(0.0, lfc_thr))
        result["UP_DE"] = uniq(ids.loc[up_mask[up_mask].index].tolist())
        result["DOWN_DE"] = uniq(ids.loc[dn_mask[dn_mask].index].tolist())
    else:
        result["UP_DE"] = []
        result["DOWN_DE"] = []

    return result

# ───────────────────────────── CLI ───────────────────────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="ORA with g:Profiler — CSV + PNG + long matrices for ALL/UP/DOWN.")
    ap.add_argument("--infile")
    ap.add_argument("--species")
    ap.add_argument("--id_column")
    ap.add_argument("--padj_thr")
    ap.add_argument("--lfc_thr")
    ap.add_argument("--sources")
    ap.add_argument("--topn")
    ap.add_argument("--outdir")
    ap.add_argument("--csv_dir")
    ap.add_argument("--fig_dir")
    args = ap.parse_args()

    padj_thr = float(args.padj_thr or 0.05)
    lfc_thr  = float(args.lfc_thr  or 0.0)
    topn     = int(args.topn or 20)

    org_code = guess_code(args.species or "Homo sapiens")
    outdir   = Path(args.outdir or ".")
    csv_dir  = Path(args.csv_dir) if args.csv_dir else (outdir / "csv")
    fig_dir  = Path(args.fig_dir) if args.fig_dir else (outdir / "figs")
    csv_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.infile, sep=None, engine="python")
    user_col = (args.id_column or "").strip() or df.columns[0]
    if user_col not in df.columns:
        raise SystemExit(f"ERROR: column '{user_col}' not found in input file.")

    id_series = series_clean_ids(df[user_col])
    if "log2FoldChange" in df.columns:
        de_lfc_map = dict(zip(id_series.astype(str), pd.to_numeric(df["log2FoldChange"], errors="coerce")))
    else:
        de_lfc_map = {}

    gene_sets = build_gene_sets(df, id_series, padj_thr=padj_thr, lfc_thr=lfc_thr)
    sources = [s.strip() for s in (args.sources or "GO:BP,GO:MF,GO:CC,KEGG").split(",") if s.strip()]
    print(f"[INFO] sources={sources}  padj_thr={padj_thr}  lfc_thr={lfc_thr}")

    for tag, genes in gene_sets.items():
        print(f"[INFO] {tag}: {len(genes)} genes")
        if not genes:
            (csv_dir / f"ORA_{tag}.csv").write_text("", encoding="utf-8")
            (csv_dir / f"ORA_{tag}_matrix_long.csv").write_text("", encoding="utf-8")
            continue

        res  = run_gprofiler(genes, org_code, sources)
        fres = filter_like_website(res, padj_thr)

        table = fres if not fres.empty else res
        if table.empty:
            (csv_dir / f"ORA_{tag}.csv").write_text("", encoding="utf-8")
        else:
            table.to_csv(csv_dir / f"ORA_{tag}.csv", index=False)

        long_df = explode_long(table if not table.empty else res, de_lfc_map)
        if long_df.empty:
            (csv_dir / f"ORA_{tag}_matrix_long.csv").write_text("", encoding="utf-8")
        else:
            long_df = long_df.drop_duplicates(subset=["term_id", "source", "gene_label"])
            long_df.to_csv(csv_dir / f"ORA_{tag}_matrix_long.csv", index=False)

        barplot(table if not table.empty else res,
                fig_dir / f"ORA_{tag}_top{topn}.png",
                f"ORA {tag} — {org_code}", topn=topn)

    print("✔ Done.")
