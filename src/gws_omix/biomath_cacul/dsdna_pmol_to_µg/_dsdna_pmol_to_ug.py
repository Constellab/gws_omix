#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, json, sys
from pathlib import Path
from typing import Dict, Any

def _float(x: str) -> float:
    if isinstance(x, str) and "," in x:
        raise argparse.ArgumentTypeError("Use '.' as decimal separator, not ','.")
    return float(x)

def dsDNA_pmol_to_ug(length_bp: float, pmol: float) -> float:
    if length_bp <= 0:
        raise ValueError("DNA Length (bp) must be > 0.")
    if pmol < 0:
        raise ValueError("DNA Amount (pmol) must be ≥ 0.")
    # µg = (pmol × N × 660 pg/pmol) × (1 µg / 1e6 pg)
    return (pmol * length_bp * 660.0) / 1e6

def _format_val_with_unit(value: float, unit: str = "µg") -> str:
    # format with up to 15 significant digits, then trim trailing zeros/dot
    s = format(value, ".15g")
    # ensure standard decimal look (no scientific unless tiny/huge)
    # (format with .15g already avoids most sci-notation for typical lab values)
    return f"{s} {unit}"

def main():
    ap = argparse.ArgumentParser(description="dsDNA: pmol → µg")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--prefix", default="dsdna_pmol_to_ug", help="Output JSON prefix")
    ap.add_argument("--length-bp", type=_float, required=True, help="DNA Length (bp)")
    ap.add_argument("--pmol", type=_float, required=True, help="DNA Amount (pmol)")
    args = ap.parse_args()

    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)

    formula = "µg = (pmol × bp × 660) / 1e6"

    try:
        val_ug = dsDNA_pmol_to_ug(args.length_bp, args.pmol)
        payload: Dict[str, Any] = {
            "tool": "dsdna_pmol_to_ug",
            "length_bp": args.length_bp,
            "pmol": args.pmol,
            "formula": formula,
            "result": _format_val_with_unit(val_ug, "µg"),
        }
    except Exception as e:
        payload = {
            "tool": "dsdna_pmol_to_ug",
            "length_bp": args.length_bp,
            "pmol": args.pmol,
            "formula": formula,
            "result": f"ERROR: {e}",
        }

    (out_dir / f"{args.prefix}.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    sys.exit(0)

if __name__ == "__main__":
    main()
