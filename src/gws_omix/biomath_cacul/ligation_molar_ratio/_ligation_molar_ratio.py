#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse, json, sys
from pathlib import Path

def _float(x: str) -> float:
    if isinstance(x, str) and "," in x:
        raise argparse.ArgumentTypeError("Use '.' as decimal separator, not ','.")
    return float(x)

def compute_insert_ng(insert_kb: float, vector_ng: float, vector_kb: float,
                      ratio_insert: float, ratio_vector: float) -> float:
    if insert_kb <= 0:
        raise ValueError("Insert Length (kb) must be > 0.")
    if vector_kb <= 0:
        raise ValueError("Vector Length (kb) must be > 0.")
    if vector_ng < 0:
        raise ValueError("Vector Amount (ng) must be ≥ 0.")
    if ratio_insert <= 0 or ratio_vector <= 0:
        raise ValueError("Ratio values (I and V) must be > 0.")
    return (ratio_insert / ratio_vector) * (insert_kb / vector_kb) * vector_ng

def main():
    ap = argparse.ArgumentParser(description="Ligations: Molar Ratio of Insert:Vector (JSON only, compact)")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--prefix", default="ligation_molar_ratio", help="Output JSON prefix")

    ap.add_argument("--insert-kb", type=_float, required=True, help="Insert Length (kb)")
    ap.add_argument("--vector-ng", type=_float, required=True, help="Vector Amount (ng)")
    ap.add_argument("--vector-kb", type=_float, required=True, help="Vector Length (kb)")
    ap.add_argument("--ratio-insert", type=_float, default=1.0, help="Insert ratio I (default 1)")
    ap.add_argument("--ratio-vector", type=_float, default=1.0, help="Vector ratio V (default 1)")
    args = ap.parse_args()

    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)
    formula = "ng_insert = (I/V) × (kb_insert / kb_vector) × ng_vector"

    try:
        val = compute_insert_ng(
            args.insert_kb, args.vector_ng, args.vector_kb,
            args.ratio_insert, args.ratio_vector
        )
        # format value with full precision but no scientific notation
        val_str = format(val, ".15g")  # e.g., 4.3010752688172
        ratio_i = int(args.ratio_insert) if args.ratio_insert.is_integer() else args.ratio_insert
        ratio_v = int(args.ratio_vector) if args.ratio_vector.is_integer() else args.ratio_vector
        result = f"{val_str} ng insert for {ratio_i}:{ratio_v} ratio"

        payload = {
            "formula": formula,
            "result": result
        }
    except Exception as e:
        payload = {
            "formula": formula,
            "result": f"ERROR: {e}"
        }

    (out_dir / f"{args.prefix}.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    sys.exit(0)

if __name__ == "__main__":
    main()
