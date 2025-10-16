#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse, json, sys
from pathlib import Path
from typing import Dict, Any

def _float(x: str) -> float:
    if isinstance(x, str) and "," in x:
        raise argparse.ArgumentTypeError("Use '.' as decimal separator, not ','.")
    return float(x)

def linear_dna_ug_to_pmol_ends(mass_ug: float, length_kb: float) -> Dict[str, float]:
    """
    Returns a dict with pmol_molecules and pmol_ends for a linear dsDNA fragment.
    """
    if length_kb <= 0:
        raise ValueError("DNA Length (kb) must be > 0.")
    if mass_ug < 0:
        raise ValueError("DNA mass (µg) must be ≥ 0.")

    length_bp = length_kb * 1000.0
    pmol_molecules = (mass_ug * 1e6) / (length_bp * 660.0)
    pmol_ends = pmol_molecules * 2.0
    return {"pmol_molecules": pmol_molecules, "pmol_ends": pmol_ends}

def main():
    ap = argparse.ArgumentParser(description="Linear DNA: µg → pmol of Ends (JSON only)")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--prefix", default="linear_dna_ug_to_pmol_ends", help="Output JSON prefix")
    ap.add_argument("--mass-ug", type=_float, required=True, help="µg of DNA")
    ap.add_argument("--length-kb", type=_float, required=True, help="DNA Length (kb)")
    args = ap.parse_args()

    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)

    payload: Dict[str, Any] = {
        "tool": "linear_dna_ug_to_pmol_ends",
        "mass_ug": args.mass_ug,
        "length_kb": args.length_kb,
        "formula": "pmol_ends = ((µg × 1e6) / (kb × 1000 × 660)) × 2"
    }

    try:
        res = linear_dna_ug_to_pmol_ends(args.mass_ug, args.length_kb)
        # ⬇⬇⬇ Modification : on n'expose plus les champs internes, seulement 'result'
        payload["result"] = f'{res["pmol_ends"]} pmol of DNA ends'
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
