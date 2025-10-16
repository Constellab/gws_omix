#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse, json, sys
from pathlib import Path
from typing import Dict

MOLARITY_FACTORS = {
    "M": 1.0,
    "MM": 1e-3,
    "UM": 1e-6,
}

VOLUME_FACTORS = {
    "L": 1.0,
    "ML": 1e-3,
    "UL": 1e-6,  
}

MASS_FACTORS_OUT = {
    "G": 1.0,
    "MG": 1e3,
    "UG": 1e6,
}

def _float(x: str) -> float:
    if isinstance(x, str) and "," in x:
        raise argparse.ArgumentTypeError("Use '.' for decimals, not ','.")
    return float(x)

def _norm_mol_unit(u: str) -> str:
    u = (u or "").strip().replace("µ", "u").upper()
    if u not in ("M", "MM", "UM"):
        raise ValueError("Concentration unit must be one of: M, mM, µM/uM.")
    return u

def _norm_vol_unit(u: str) -> str:
    u = (u or "").strip().replace("µ", "u").upper()
    if u not in ("L", "ML", "UL"):
        raise ValueError("Volume unit must be one of: L, mL, µL/uL.")
    return u

def _norm_mass_unit_out(u: str | None) -> str:
    if not u:
        return "G"
    u = (u or "").strip().replace("µ", "u").upper()
    if u not in ("G", "MG", "UG"):
        raise ValueError("Output mass unit must be one of: g, mg, µg/uG.")
    return u

def _unit_lower_mass(u_norm: str) -> str:
    return {"G": "g", "MG": "mg", "UG": "ug"}[u_norm]

def main():
    ap = argparse.ArgumentParser(description="Molarity Calculator (JSON only): mass = (M × V) × MW.")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--prefix", default="molarity", help="Output JSON prefix")

    ap.add_argument("--mw", type=_float, required=True, help="Molecular weight (g/mol)")
    ap.add_argument("--molarity", type=_float, required=True, help="Final (target) molarity (numeric)")
    ap.add_argument("--mol-unit", required=True, help="Unit for molarity: M, mM, µM/uM")
    ap.add_argument("--final-vol", type=_float, required=True, help="Final volume (numeric)")
    ap.add_argument("--vol-unit", required=True, help="Unit for final volume: L, mL, µL/uL")

    ap.add_argument("--mass-unit-out", required=False, help="Output mass unit: g, mg, µg/uG (default: g)")

    args = ap.parse_args()
    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)

    payload: Dict = {
        "tool": "molarity_calculator",
        "inputs": {
            "molecular_weight_g_per_mol": args.mw,
            "final_concentration": args.molarity,
            "final_concentration_unit": args.mol_unit,
            "final_volume": args.final_vol,
            "final_volume_unit": args.vol_unit,
            "mass_unit_out": args.mass_unit_out or "g"
        },
        "formula": "moles = M × V ; mass(g) = moles × MW",
    }

    try:
        if args.mw <= 0:
            raise ValueError("Molecular weight must be > 0.")
        if args.molarity <= 0:
            raise ValueError("Final molarity must be > 0.")
        if args.final_vol <= 0:
            raise ValueError("Final volume must be > 0.")

        uM = _norm_mol_unit(args.mol_unit)
        uV = _norm_vol_unit(args.vol_unit)
        uMassOut = _norm_mass_unit_out(args.mass_unit_out)

        M_in_M = args.molarity * MOLARITY_FACTORS[uM]  # mol/L
        V_in_L = args.final_vol * VOLUME_FACTORS[uV]   # L

        moles = M_in_M * V_in_L                        # mol
        mass_g = moles * args.mw                       # g

        # Convert to desired output mass unit (display only)
        mass_out = mass_g * MASS_FACTORS_OUT[uMassOut]
        unit_out = _unit_lower_mass(uMassOut)          # "g" | "mg" | "ug"

        payload["result"] = {
            "results": f"{mass_out} {unit_out}",
            "message": f"weigh {mass_out} {unit_out} of compound"
        }

    except Exception as e:
        payload["error"] = str(e)

    (out_dir / f"{args.prefix}.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    # Always rc=0 to surface JSON even on input mistakes
    sys.exit(0)

if __name__ == "__main__":
    main()
