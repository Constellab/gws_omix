#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse, json, sys
from pathlib import Path
from typing import Dict, Any

def _float(x: str) -> float:
    if isinstance(x, str) and "," in x:
        raise argparse.ArgumentTypeError("Utiliser '.' comme séparateur décimal, pas ','.")
    return float(x)

def dsDNA_ug_to_pmol(length_bp: float, mass_ug: float) -> float:
    if length_bp <= 0:
        raise ValueError("DNA Length (bp) doit être > 0.")
    if mass_ug < 0:
        raise ValueError("DNA Amount (µg) doit être ≥ 0.")
    return (mass_ug * 1e6) / (length_bp * 660.0)

def main():
    ap = argparse.ArgumentParser(description="dsDNA: µg → pmol (JSON only)")
    ap.add_argument("--out", required=True, help="Dossier de sortie")
    ap.add_argument("--prefix", default="dsdna_ug_to_pmol", help="Préfixe fichier JSON")
    ap.add_argument("--length-bp", type=_float, required=True, help="DNA Length (bp)")
    ap.add_argument("--mass-ug", type=_float, required=True, help="DNA Amount (µg)")
    args = ap.parse_args()

    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)
    payload: Dict[str, Any] = {
        "tool": "dsdna_ug_to_pmol",
        "length_bp": args.length_bp,
        "mass_ug": args.mass_ug,
        "formula": "pmol = (µg × 1e6) / (bp × 660)"
    }

    try:
        # ⬇⬇⬇ CHANGEMENT ICI UNIQUEMENT : on écrit 'result' "<valeur> pmol"
        val = dsDNA_ug_to_pmol(args.length_bp, args.mass_ug)
        payload["result"] = f"{val} pmol"
    except Exception as e:
        payload["error"] = str(e)
        # on écrit quand même le JSON d'erreur
        (out_dir / f"{args.prefix}.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
        print(json.dumps(payload, indent=2, ensure_ascii=False))
        sys.exit(2)

    # JSON only
    (out_dir / f"{args.prefix}.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    sys.exit(0)

if __name__ == "__main__":
    main()
