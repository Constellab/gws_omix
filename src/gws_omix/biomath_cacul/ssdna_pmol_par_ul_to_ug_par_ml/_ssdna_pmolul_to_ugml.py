#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse, json, sys
from pathlib import Path
from typing import Dict, Any

def _float(x: str) -> float:
    if isinstance(x, str) and "," in x:
        raise argparse.ArgumentTypeError("Use '.' as decimal separator, not ','.")
    return float(x)

def ssdna_pmolul_to_ugml(length_nt: float, conc_pmolul: float) -> float:
    if length_nt <= 0:
        raise ValueError("Oligo Length (bases) must be > 0.")
    if conc_pmolul < 0:
        raise ValueError("Concentration (pmol/µl) must be ≥ 0.")
    # µg/ml = (pmol/µl × N × 330) / 1000
    return (conc_pmolul * length_nt * 330.0) / 1000.0

def main():
    ap = argparse.ArgumentParser(description="ssDNA: pmol/µl → µg/ml (JSON only)")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--prefix", default="ssdna_pmolul_to_ugml", help="Output JSON prefix")
    ap.add_argument("--length-nt", type=_float, required=True, help="Oligo Length (bases)")
    ap.add_argument("--conc-pmolul", type=_float, required=True, help="Concentration (pmol/µl)")
    args = ap.parse_args()

    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)

    payload: Dict[str, Any] = {
        "tool": "ssdna_pmolul_to_ugml",
        "length_nt": args.length_nt,
        "conc_pmolul": args.conc_pmolul,
        "formula": "µg/ml = (pmol/µl × N × 330) / 1000"
    }

    try:
        val = ssdna_pmolul_to_ugml(args.length_nt, args.conc_pmolul)
        # ✅ changement: un seul champ 'result' avec l'unité
        payload["result"] = f"{val} µg/ml"
    except Exception as e:
        payload["error"] = str(e)
        (out_dir / f"{args.prefix}.json").write_text(
            json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        print(json.dumps(payload, indent=2, ensure_ascii=False))
        sys.exit(2)

    (out_dir / f"{args.prefix}.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    sys.exit(0)

if __name__ == "__main__":
    main()
