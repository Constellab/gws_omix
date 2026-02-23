#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import re
import time
import random
import urllib.request
import urllib.error
from typing import Dict, List, Tuple

import pandas as pd


def detect_sep(path: str) -> str:
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        sample = fh.read(65536)
    candidates = [",", "\t", ";", "|"]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="".join(candidates))
        return dialect.delimiter
    except Exception:
        header = sample.splitlines()[0] if sample else ""
        counts = {d: header.count(d) for d in candidates}
        return max(counts, key=counts.get) if counts else ","


def read_table(path: str) -> pd.DataFrame:
    sep = detect_sep(path)
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, engine="python")


def clean_entrez(x: str) -> str:
    """
    Accept:
      - "318213.0", " 318213 ", "GeneID:318213"
    Return:
      - "318213" or "" if invalid
    """
    x = str(x).strip()
    x = re.sub(r"\.0$", "", x)
    x = re.sub(r"[^\d]", "", x)
    return x if re.fullmatch(r"\d+", x) else ""


def chunk_entries(entries: List[str], max_items: int, max_url_len: int) -> List[List[str]]:
    chunks: List[List[str]] = []
    cur: List[str] = []
    cur_len = 0
    for e in entries:
        e = str(e).strip()
        if not e:
            continue
        add = len(e) + (1 if cur else 0)
        if cur and (len(cur) >= max_items or (cur_len + add) > max_url_len):
            chunks.append(cur)
            cur = [e]
            cur_len = len(e)
        else:
            cur.append(e)
            cur_len += add
    if cur:
        chunks.append(cur)
    return chunks


def fetch_text(url: str, timeout: int, max_retries: int) -> str:
    attempt = 1
    while True:
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and attempt <= max_retries:
                sleep_s = min((2 ** attempt) * 0.25, 20.0) + random.uniform(0, 0.3)
                time.sleep(sleep_s)
                attempt += 1
                continue
            raise
        except Exception:
            if attempt <= max_retries:
                sleep_s = min((2 ** attempt) * 0.25, 20.0) + random.uniform(0, 0.3)
                time.sleep(sleep_s)
                attempt += 1
                continue
            raise


def _split_line_2cols(line: str) -> Tuple[str, str]:
    line = line.strip()
    if "\t" in line:
        a, b = line.split("\t", 1)
        return a.strip(), b.strip()
    parts = line.split(None, 1)
    if len(parts) < 2:
        return "", ""
    return parts[0].strip(), parts[1].strip()


def kegg_conv_entrez_to_species(
    specie: str,
    entrez_ids: List[str],
    max_items_per_query: int,
    max_url_len: int,
    timeout_sec: int,
    max_retries: int,
    sleep_between_calls_sec: float,
) -> Dict[str, str]:
    """
    GET https://rest.kegg.jp/conv/<org>/ncbi-geneid:<id>+...
    Returns mapping: entrez_id -> "<org>:<gene>"
    """
    ids = [str(x).strip() for x in entrez_ids if str(x).strip()]
    ids = sorted(set(ids))
    if not ids:
        return {}

    prefixed = [f"ncbi-geneid:{x}" for x in ids]
    chunks = chunk_entries(prefixed, max_items_per_query, max_url_len)

    mapping: Dict[str, str] = {}
    for sub in chunks:
        url = f"https://rest.kegg.jp/conv/{specie}/" + "+".join(sub)
        raw = fetch_text(url, timeout=timeout_sec, max_retries=max_retries)

        for line in raw.splitlines():
            if not line.strip():
                continue
            a, b = _split_line_2cols(line)
            if not a or not b:
                continue

            if a.startswith(specie + ":") and b.startswith("ncbi-geneid:"):
                eid = b.split(":", 1)[1]
                mapping[eid] = a
            elif b.startswith(specie + ":") and a.startswith("ncbi-geneid:"):
                eid = a.split(":", 1)[1]
                mapping[eid] = b

        time.sleep(sleep_between_calls_sec + random.uniform(0, 0.05))

    return mapping


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deg", required=True)
    ap.add_argument("--specie", required=True)
    ap.add_argument("--col_entrez", required=True)
    ap.add_argument(
        "--foldchange_cols",
        default="",
        help="Comma-separated list of fold-change columns, e.g. 'FC1,FC2,log2FoldChange'",
    )
    ap.add_argument("--min_genes_required", type=int, default=1)
    ap.add_argument("--out_gene_kegg", required=True)
    ap.add_argument("--out_used", required=True)
    args = ap.parse_args()

    TIMEOUT_SEC = 180
    MAX_RETRIES = 6
    SLEEP_BETWEEN_CALLS_SEC = 0.15
    MAX_ITEMS_PER_QUERY = 80
    MAX_URL_LEN = 1800

    df = read_table(args.deg)
    if df.empty:
        raise SystemExit("Input DEG is empty.")

    if args.col_entrez not in df.columns:
        raise SystemExit(f"Entrez column '{args.col_entrez}' not found in input columns: {list(df.columns)}")

    # ---- Clean & validate Entrez column ----
    raw = df[args.col_entrez].astype(str).str.strip()
    cleaned = raw.map(clean_entrez)
    n_valid = int((cleaned != "").sum())

    # Hard check: if user provided non-numeric IDs, stop with a clear message
    if n_valid == 0:
        sample_bad = raw[raw != ""].head(10).tolist()
        raise SystemExit(
            "Selected column does not contain valid NCBI Entrez Gene IDs (numbers). "
            "Please run Gene ID conversion (OmiX) to ENTREZGENE and re-run with that numeric column. "
            f"Example values seen: {sample_bad}"
        )

    # Keep only rows with valid Entrez
    df["_entrez"] = cleaned
    df = df[df["_entrez"] != ""].copy()

    # ---- Convert Entrez -> KEGG gene ----
    unique_ids = sorted(set(df["_entrez"].tolist()))
    conv = kegg_conv_entrez_to_species(
        args.specie,
        unique_ids,
        MAX_ITEMS_PER_QUERY,
        MAX_URL_LEN,
        TIMEOUT_SEC,
        MAX_RETRIES,
        SLEEP_BETWEEN_CALLS_SEC,
    )

    df["kegg_full"] = df["_entrez"].map(lambda x: conv.get(x, ""))
    df["kegg_short"] = df["kegg_full"].map(lambda x: x.split(":", 1)[1] if ":" in x else "")

    mapped = df[df["kegg_short"] != ""].copy()
    if mapped.shape[0] < args.min_genes_required:
        raise SystemExit(f"Too few genes mapped to KEGG ({mapped.shape[0]} < {args.min_genes_required}).")

    # ---- Fold-change columns (optional) ----
    fc_cols = [c.strip() for c in (args.foldchange_cols or "").split(",") if c.strip()]
    for c in fc_cols:
        if c not in mapped.columns:
            raise SystemExit(f"Fold-change column '{c}' not found in input.")

    out = mapped[["kegg_short"] + fc_cols].copy()

    # Make FC numeric where possible (Pathview needs numeric)
    for c in fc_cols:
        out[c] = pd.to_numeric(out[c], errors="coerce")

    out.to_csv(args.out_gene_kegg, index=False)

    # Debug / trace file (not an output of the task)
    df.to_csv(args.out_used, index=False)

    print(f"Valid Entrez IDs in input: {n_valid}")
    print(f"Mapped rows to KEGG genes: {mapped.shape[0]}/{df.shape[0]}")


if __name__ == "__main__":
    main()
