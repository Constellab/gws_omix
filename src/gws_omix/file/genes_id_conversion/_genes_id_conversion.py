#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
g:Profiler ID conversion CLI

What it does
------------
1) Reads a CSV/TSV (auto-detects delimiter unless --sep is passed).
2) Calls g:Profiler gconvert (POST JSON) to map identifiers.
3) Writes:
   - <out_prefix>_converted.csv  : the tidy conversion table returned by g:Profiler
   - <out_prefix>_annotated.csv/tsv : your original table + new columns (converted, name, description)
   - <out_prefix>_summary.json  : a small JSON summary printed to stdout too

Notes
-----
- The column used for conversion is chosen with --id_column (default: first column).
- --organism accepts either a g:Profiler code (e.g., "hsapiens") or a scientific name (e.g., "Homo sapiens").
- --numeric_ns can be passed to treat bare numeric IDs correctly (e.g., ENTREZGENE_ACC). Omit or use empty to disable.
"""

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
import requests

API_BASE = "https://biit.cs.ut.ee/gprofiler/api/convert/convert/"  # POST JSON


def _detect_sep(fp: str, sep_arg: Optional[str]) -> str:
    """
    Detect the delimiter of a text table, unless the user explicitly provided --sep.
    Returns a single-character separator. Defaults to comma.
    """
    if sep_arg is not None and len(sep_arg) > 0:
        return sep_arg

    with open(fp, "r", encoding="utf-8", errors="replace") as fh:
        sample = fh.read(64_000)

    candidates = [",", "\t", ";", "|"]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="".join(candidates))
        return dialect.delimiter
    except Exception:
        header = sample.splitlines()[0] if sample else ""
        counts = {d: header.count(d) for d in candidates}
        return max(counts, key=counts.get) if counts else ","


def _read_table(fp: str, sep: Optional[str]) -> pd.DataFrame:
    """
    Robust reader:
      - If sep provided, use it.
      - Else sniff delimiter from a sample; fall back to best of [',','\\t',';','|'].
      - Force dtype=str (avoid numeric coercion / gene symbol mangling).
    """
    if sep is not None:
        return pd.read_csv(fp, sep=sep, dtype=str, keep_default_na=False)

    with open(fp, "r", encoding="utf-8", errors="replace") as fh:
        sample = fh.read(64_000)

    candidates = [",", "\t", ";", "|"]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="".join(candidates))
        guessed = dialect.delimiter
    except Exception:
        header = sample.splitlines()[0] if sample else ""
        counts = {d: header.count(d) for d in candidates}
        guessed = max(counts, key=counts.get) if counts else ","

    return pd.read_csv(
        fp,
        sep=guessed,
        dtype=str,
        keep_default_na=False,
        engine="python",
        quoting=csv.QUOTE_MINIMAL,
    )


def _to_org_code(s: str) -> str:
    """
    Accept either a g:Profiler code ('hsapiens') or a scientific name ('Homo sapiens').
    'Genus species' -> 'gspecies' lowercase; otherwise return lowercase as-is.
    """
    val = (s or "").strip()
    if not val:
        raise SystemExit("organism is required")
    if " " in val or "_" in val:
        parts = re.split(r"[\s_]+", val)
        if len(parts) < 2:
            raise SystemExit(f"Unrecognized organism string: {s}")
        return (parts[0][0] + parts[1]).lower()
    return val.lower()


def _gconvert_post(org_code: str, ids: List[str], target_ns: str, numeric_ns: Optional[str]) -> List[Dict[str, Any]]:
    """
    Call g:Profiler gconvert via POST JSON.
    - Required: organism, query[], target
    - Optional: numeric_ns (e.g., 'ENTREZGENE_ACC') to treat bare numeric IDs properly.
    Returns a list of mapping rows (JSON objects).
    """
    payload: Dict[str, Any] = {
        "organism": org_code,
        "query": ids,        # send list, not a single string
        "target": target_ns, # g:Profiler calls this field 'target'
    }
    if numeric_ns:
        payload["numeric_ns"] = numeric_ns

    resp = requests.post(API_BASE, json=payload, timeout=60)
    if resp.status_code != 200:
        raise RuntimeError(f"g:Profiler gconvert API error {resp.status_code}: {resp.text}")

    data = resp.json()
    # Response shape can be either a dict with 'result' or a raw list; handle both.
    rows = data.get("result") if isinstance(data, dict) else data
    if rows is None:
        rows = []
    if not isinstance(rows, list):
        raise RuntimeError("Unexpected response from g:Profiler gconvert API.")
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description="g:Profiler ID conversion")
    ap.add_argument("--infile", required=True)
    ap.add_argument("--id_column", default=None,
                    help="Optional: column name to convert (defaults to first column).")
    ap.add_argument("--organism", required=True)
    ap.add_argument("--target_ns", required=True)
    ap.add_argument("--numeric_ns", default="",
                    help="Optional: how to treat *bare numeric* IDs (e.g., ENTREZGENE_ACC).")
    ap.add_argument("--out_prefix", default="ID_CONVERT")
    ap.add_argument("--sep", default=None,
                    help="Optional explicit separator (e.g., ',' or '\\t').")
    args = ap.parse_args()

    # Detect and remember the input separator to reuse it for the annotated file
    input_sep = _detect_sep(args.infile, args.sep)

    # 1) Load input
    df = _read_table(args.infile, input_sep)
    if df is None or df.empty:
        raise SystemExit("Input is empty.")

    # 2) Choose ID column (optional → falls back to first column)
    idcol = args.id_column or df.columns[0]
    if idcol not in df.columns:
        raise SystemExit(f"Column '{idcol}' not found. Available: {list(df.columns)}")

    # 3) Extract IDs as clean strings
    ids = (
        df[idcol]
        .astype(str)
        .str.strip()
        .replace({"": pd.NA})
        .dropna()
        .tolist()
    )
    if not ids:
        raise SystemExit("No IDs to convert after dropping empty values.")

    # 4) Convert via HTTP API (supports numeric_ns)
    org_code = _to_org_code(args.organism)
    numeric_ns = (args.numeric_ns or "").strip() or None

    rows = _gconvert_post(org_code, ids, args.target_ns, numeric_ns)

    # 5) Normalize columns to a stable set
    # Keep tidy names: incoming, converted, name, description
    norm = []
    for r in rows:
        incoming = r.get("incoming", r.get("input"))
        converted = r.get("converted", r.get("name"))  # 'name' is often the mapped ID
        # keep an extra alias if available
        name = r.get("name") if "converted" in r else r.get("converted")
        description = r.get("description")
        norm.append({
            "incoming": incoming,
            "converted": converted,
            "name": name,
            "description": description,
        })

    mapped = pd.DataFrame(norm, columns=["incoming", "converted", "name", "description"])

    # 6) Write converted-only table (always CSV)
    out_prefix = Path(args.out_prefix)
    converted_fp = f"{out_prefix}_converted.csv"
    mapped.to_csv(converted_fp, index=False)

    # 7) Build annotated table = original DF + conversion columns
    #    - deduplicate by incoming to avoid exploding merges (first match wins)
    #    - merge on the original id column
    if not mapped.empty:
        dedup = (
            mapped.drop_duplicates(subset=["incoming"], keep="first")
                  .rename(columns={"incoming": idcol})
        )
        # Avoid duplicating the join key if present
        drop_cols = [c for c in [idcol] if c in dedup.columns and c != idcol]
        annotated = df.merge(dedup.drop(columns=drop_cols, errors="ignore"), how="left", on=idcol)
    else:
        # No mappings returned; just copy original and append empty columns
        annotated = df.copy()
        for col in ["converted", "name", "description"]:
            if col not in annotated.columns:
                annotated[col] = pd.NA

    # 8) Save annotated file using the same delimiter as the input
    ext = "tsv" if input_sep == "\t" else "csv"
    annotated_fp = f"{out_prefix}_annotated.{ext}"
    annotated.to_csv(annotated_fp, index=False, sep=input_sep)

    # 9) Summary JSON (useful for logs)
    n_mapped_effective = int(mapped["converted"].notna().sum()) if not mapped.empty else 0
    summary = {
        "n_input": len(ids),
        "n_mapped": n_mapped_effective,
        "organism": args.organism,
        "organism_code": org_code,
        "target_namespace": args.target_ns,
        "numeric_ns": numeric_ns or "",
        "input_file": str(args.infile),
        "id_column": idcol,
        "converted_table": converted_fp,
        "annotated_table": annotated_fp,
    }
    Path(f"{out_prefix}_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
