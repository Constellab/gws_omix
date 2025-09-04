#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, csv, json, re
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
import requests


API_BASE = "https://biit.cs.ut.ee/gprofiler/api/convert/convert/"  # POST JSON


def _read_table(fp: str, sep: str | None) -> pd.DataFrame:
    """
    Robust reader:
      - If --sep provided, use it.
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


def _gconvert_post(org_code: str, ids: List[str], target_ns: str, numeric_ns: str | None) -> List[Dict[str, Any]]:
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


def main():
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

    # 1) Load input
    df = _read_table(args.infile, args.sep)
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
    # g:Profiler typically returns fields like: input / incoming, name (mapped ID), description, ...
    # We keep your preferred tidy names: incoming, converted, name, description
    norm = []
    for r in rows:
        incoming = r.get("incoming", r.get("input"))
        converted = r.get("converted", r.get("name"))  # 'name' is often the mapped ID
        name = r.get("name") if "converted" in r else r.get("converted")  # keep an extra alias if available
        description = r.get("description")
        norm.append({
            "incoming": incoming,
            "converted": converted,
            "name": name,
            "description": description,
        })

    mapped = pd.DataFrame(norm, columns=["incoming", "converted", "name", "description"])

    # 6) Write outputs
    out_prefix = Path(args.out_prefix)
    mapped.to_csv(f"{out_prefix}_converted.csv", index=False)

    # 7) Summary JSON (useful for logs)
    # Success count = non-empty 'converted'
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
    }
    Path(f"{out_prefix}_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
