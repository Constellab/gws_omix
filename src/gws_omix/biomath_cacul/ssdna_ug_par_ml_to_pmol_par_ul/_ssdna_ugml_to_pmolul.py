#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse, json, sys
from pathlib import Path
from typing import Dict, Any

def _float(x: str) -> float:
    if isinstance(x, str) and "," in x:
        raise argparse.ArgumentTypeError("Use '.' as decimal separator, not ','.")
    return float(x)

def ssdna_ugml_to_pmolul(length_nt: float, conc_ugml: float) -> float:
    if length_nt <= 0:
        raise ValueError("Oligo Length (bases) must be > 0.")
    if conc_ugml < 0:
        raise ValueError("Concentration (µg/ml) must be ≥ 0.")
    # pmol/µl = (µg/ml × 1000) / (N × 330)
    return (conc_ugml * 1000.0) / (length_nt * 330.0)

def main():
    ap = argparse.ArgumentParser(description="ssDNA: µg/ml → pmol/µl (JSON only)")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--prefix", default="ssdna_ugml_to_pmolul", help="Output JSON prefix")
    ap.add_argument("--length-nt", type=_float, required=True, help="Oligo Length (bases)")
    ap.add_argument("--conc-ugml", type=_float, required=True, help="Concentration (µg/ml)")
    args = ap.parse_args()

    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)

    payload: Dict[str, Any] = {
        "tool": "ssdna_ugml_to_pmolul",
        "length_nt": args.length_nt,
        "conc_ugml": args.conc_ugml,
        "formula": "pmol/µl = (µg/ml × 1000) / (N × 330)"
    }

    try:
        val = ssdna_ugml_to_pmolul(args.length_nt, args.conc_ugml)
        # ✅ changement: un seul champ 'result' avec l'unité
        payload["result"] = f"{val} pmol/µl"
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
