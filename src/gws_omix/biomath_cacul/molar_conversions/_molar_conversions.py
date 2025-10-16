#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse, json, sys
from pathlib import Path
from typing import Dict, Any, Tuple

CHOICES = ["Protein Size (kDa)", "pmol of Protein", "μg of Protein"]

def _float(x: str) -> float:
    if isinstance(x, str) and "," in x:
        raise argparse.ArgumentTypeError("Use '.' for decimals, not ','.")
    return float(x)

def normalize_mode(label: str) -> str:
    """
    Map UI label to internal token: 'kda', 'pmol', or 'ug'.
    """
    s = (label or "").strip().lower()
    if s == "protein size (kda)":
        return "kda"
    if s == "pmol of protein":
        return "pmol"
    if s in ("μg of protein", "ug of protein"):
        return "ug"
    # soft aliases
    if s in ("size", "kda"): return "kda"
    if s in ("pmol",): return "pmol"
    if s in ("ug", "µg"): return "ug"
    raise ValueError('calculate-for must be one of: "Protein Size (kDa)", "pmol of Protein", "μg of Protein".')

def compute(mode: str, size_kda: float=None, pmol: float=None, mass_ug: float=None) -> Tuple[float, str, str]:
    """
    Compute the requested quantity and return (value, unit, formula_text).
      mode = 'ug'  : need kDa and pmol          → µg = (kDa × pmol) / 1000
      mode = 'pmol': need kDa and µg            → pmol = (µg × 1000) / kDa
      mode = 'kda' : need µg and pmol           → kDa = (µg × 1000) / pmol
    """
    if mode == "ug":
        if size_kda is None or pmol is None:
            raise ValueError("To calculate µg of Protein, provide Protein Size (kDa) and pmol of Protein.")
        if size_kda <= 0 or pmol < 0:
            raise ValueError("Protein Size (kDa) must be > 0 and pmol must be ≥ 0.")
        val = (size_kda * pmol) / 1000.0
        return val, "µg", "µg = (kDa × pmol) / 1000"

    if mode == "pmol":
        if size_kda is None or mass_ug is None:
            raise ValueError("To calculate pmol of Protein, provide Protein Size (kDa) and µg of Protein.")
        if size_kda <= 0 or mass_ug < 0:
            raise ValueError("Protein Size (kDa) must be > 0 and µg must be ≥ 0.")
        val = (mass_ug * 1000.0) / size_kda
        return val, "pmol", "pmol = (µg × 1000) / kDa"

    if mode == "kda":
        if pmol is None or mass_ug is None:
            raise ValueError("To calculate Protein Size (kDa), provide µg of Protein and pmol of Protein.")
        if pmol <= 0 or mass_ug < 0:
            raise ValueError("pmol must be > 0 and µg must be ≥ 0.")
        val = (mass_ug * 1000.0) / pmol
        return val, "kDa", "kDa = (µg × 1000) / pmol"

    raise ValueError("Unknown mode.")

def main():
    ap = argparse.ArgumentParser(description="Molar Conversions (Protein) — JSON only")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--prefix", default="molar_conversions", help="Output JSON prefix")
    ap.add_argument("--calculate-for", required=True, help='One of: "Protein Size (kDa)" | "pmol of Protein" | "μg of Protein"')
    ap.add_argument("--size-kda", type=_float, help="Protein Size (kDa)")
    ap.add_argument("--pmol", type=_float, help="pmol of Protein")
    ap.add_argument("--mass-ug", type=_float, help="µg of Protein")
    args = ap.parse_args()

    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)

    payload: Dict[str, Any] = {
        "tool": "molar_conversions",
        "calculate_for": args.calculate_for,
        "size_kda": args.size_kda,
        "pmol": args.pmol,
        "mass_ug": args.mass_ug
    }

    try:
        mode = normalize_mode(args.calculate_for)
        value, unit, formula = compute(mode, size_kda=args.size_kda, pmol=args.pmol, mass_ug=args.mass_ug)
        payload["formula"] = formula
        # 6 sig figs, no scientific unless needed
        payload["result"] = f"{format(value, '.6g')} {unit}"
    except Exception as e:
        payload["error"] = str(e)
        payload["result"] = f"ERROR: {e}"

    (out_dir / f"{args.prefix}.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    sys.exit(0)

if __name__ == "__main__":
    main()
