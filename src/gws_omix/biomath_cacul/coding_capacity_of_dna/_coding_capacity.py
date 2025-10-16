#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse, json, sys
from pathlib import Path
from typing import Dict, Any

AA_TO_KDA = 0.11

def _float(x: str) -> float:
    if isinstance(x, str) and "," in x:
        raise argparse.ArgumentTypeError("Use '.' as decimal separator, not ','.")
    return float(x)

def count_non_none(*vals) -> int:
    return sum(v is not None for v in vals)

def compute_from_bp(dna_bp: float) -> Dict[str, float]:
    aa = dna_bp / 3.0
    kda = aa * AA_TO_KDA
    return {"DNA Length (bp)": dna_bp, "Protein Length (amino acids)": aa, "Protein Size (kDa)": kda}

def compute_from_aa(aa: float) -> Dict[str, float]:
    dna_bp = 3.0 * aa
    kda = aa * AA_TO_KDA
    return {"DNA Length (bp)": dna_bp, "Protein Length (amino acids)": aa, "Protein Size (kDa)": kda}

def compute_from_kda(kda: float) -> Dict[str, float]:
    aa = kda / AA_TO_KDA
    dna_bp = 3.0 * aa
    return {"DNA Length (bp)": dna_bp, "Protein Length (amino acids)": aa, "Protein Size (kDa)": kda}

def main():
    ap = argparse.ArgumentParser(description="Coding Capacity of DNA — JSON only (enter exactly one parameter).")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--prefix", default="coding_capacity", help="Output JSON prefix")

    # Exactly one of these must be provided
    ap.add_argument("--dna-bp", type=_float, help="DNA Length (bp)")
    ap.add_argument("--aa", type=_float, help="Protein Length (amino acids)")
    ap.add_argument("--protein-kda", type=_float, help="Protein Size (kDa)")

    args = ap.parse_args()
    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)

    provided = count_non_none(args.dna_bp, args.aa, args.protein_kda)

    # Build input section using pretty keys (even if missing)
    input_block = {
        "DNA Length (bp)": args.dna_bp,
        "Protein Length (amino acids)": args.aa,
        "Protein Size (kDa)": args.protein_kda,
    }

    payload: Dict[str, Any] = {
        "tool": "coding_capacity",
        "input": input_block,
        "formulas": {
            "Protein Length (amino acids) from DNA Length (bp)": "aa = bp / 3",
            "DNA Length (bp) from Protein Length (amino acids)": "bp = 3 × aa",
            "Protein Size (kDa) from Protein Length (amino acids)": "kDa = aa × 0.11",
        },
        "note": "Provide exactly one input (DNA Length (bp) OR Protein Length (amino acids) OR Protein Size (kDa))."
    }

    try:
        if provided == 0:
            raise ValueError("Provide exactly one parameter: --dna-bp OR --aa OR --protein-kda.")
        if provided > 1:
            raise ValueError("Only one parameter must be provided at a time (not multiple).")

        if args.dna_bp is not None:
            if args.dna_bp <= 0:
                raise ValueError("DNA Length (bp) must be > 0.")
            res = compute_from_bp(args.dna_bp)
            path = "from DNA Length (bp)"
        elif args.aa is not None:
            if args.aa <= 0:
                raise ValueError("Protein Length (amino acids) must be > 0.")
            res = compute_from_aa(args.aa)
            path = "from Protein Length (amino acids)"
        else:
            if args.protein_kda <= 0:
                raise ValueError("Protein Size (kDa) must be > 0.")
            res = compute_from_kda(args.protein_kda)
            path = "from Protein Size (kDa)"

        payload["path"] = path
        payload["results"] = res

    except Exception as e:
        payload["error"] = str(e)

    (out_dir / f"{args.prefix}.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    # Always 0 so the JSON is surfaced even on input mistakes
    sys.exit(0)

if __name__ == "__main__":
    main()
