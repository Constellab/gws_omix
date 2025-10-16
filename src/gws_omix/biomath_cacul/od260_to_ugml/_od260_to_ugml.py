#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse, json, sys
from pathlib import Path

FACTORS = {
    "DNA": 50.0,
    "RNA": 40.0,
    "Single-stranded-DNA": 35.0,
    "Single-stranded-Oligo": 20.0
}

def _float(x: str) -> float:
    try:
        return float(x.replace(",", "."))
    except Exception:
        raise argparse.ArgumentTypeError("Use '.' for decimals.")

def normalize_sample_type(s: str) -> str:
    """Normalise le type d'échantillon avec tirets (et accepte abréviations)."""
    s = (s or "").strip()
    low = s.lower()

    if low in ("dna", "dsdna"):
        return "DNA"
    if low in ("rna", "ssrna"):
        return "RNA"
    if low in ("ssdna", "single-stranded-dna", "single stranded dna"):
        return "Single-stranded-DNA"
    if low in ("ssoligo", "single-stranded-oligo", "single stranded oligo", "oligo"):
        return "Single-stranded-Oligo"
    if low == "single-stranded":
        # Valeur tronquée — on ne laisse plus d'erreur bloquante
        # on suppose DNA comme fallback
        return "Single-stranded-DNA"

    raise ValueError(
        "Sample Type must be one of: DNA, RNA, Single-stranded-DNA, Single-stranded-Oligo."
    )

def compute_ugml(od260: float, sample_type: str) -> float:
    return od260 * FACTORS[sample_type]

def main():
    parser = argparse.ArgumentParser(description="OD260 → µg/ml (Promega version tirets)")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--prefix", default="od260_to_ugml", help="Output JSON prefix")
    parser.add_argument("--od260", type=_float, required=True, help="OD260 value")
    parser.add_argument("--sample-type", required=True, help="DNA | RNA | Single-stranded-DNA | Single-stranded-Oligo")

    args = parser.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    formula = "µg/ml = OD260 × conversion factor"

    try:
        s_type = normalize_sample_type(args.sample_type)
        result = compute_ugml(args.od260, s_type)
        payload = {
            "tool": "od260_to_ugml",
            "od260": args.od260,
            "sample_type": s_type,
            "factor_ugml_per_od": FACTORS[s_type],
            "formula": formula,
            "result": f"{result:.2f} µg/ml"
        }
    except Exception as e:
        payload = {
            "tool": "od260_to_ugml",
            "od260": args.od260,
            "sample_type": args.sample_type,
            "formula": formula,
            "result": f"ERROR: {str(e)}"
        }

    # Écriture JSON
    output = out_dir / f"{args.prefix}.json"
    output.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    sys.exit(0)

if __name__ == "__main__":
    main()
