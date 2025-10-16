#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Dilution Calculator — JSON worker
---------------------------------

Calcule le volume de solution mère (V1) pour préparer une solution finale (C2, V2) :
    C1 · V1 = C2 · V2  ⇒  V1 = (C2 · V2) / C1

Unités supportées :
  - Concentrations : M, mM, µM/uM, nM, pM
  - Volumes        : L, mL, µL/uL

Nouveau :
  - --stock-vol-unit-out : unité d'affichage de V1 (L, mL, µL/uL). Si omis, V1 est affiché
    dans la même unité que le volume final (comportement historique).

Sortie (section 'result' uniquement) :
  - results : "<valeur> <unité_minuscule>"
  - message : "add <valeur> <unité_minuscule> of stock solution"
"""

import argparse, json, sys
from pathlib import Path
from typing import Dict

CONC_FACTORS = {"M": 1.0, "MM": 1e-3, "UM": 1e-6, "NM": 1e-9, "PM": 1e-12}
VOL_FACTORS  = {"L": 1.0, "ML": 1e-3, "UL": 1e-6}

def _float(x: str) -> float:
    if isinstance(x, str) and "," in x:
        raise argparse.ArgumentTypeError("Use '.' for decimals, not ','.")
    return float(x)

def _norm_conc_unit(u: str) -> str:
    u = (u or "").strip().replace("µ", "u").upper()
    if u not in ("M", "MM", "UM", "NM", "PM"):
        raise ValueError("Concentration unit must be one of: M, mM, µM/uM, nM, pM.")
    return u

def _norm_vol_unit(u: str) -> str:
    u = (u or "").strip().replace("µ", "u").upper()
    if u not in ("L", "ML", "UL"):
        raise ValueError("Volume unit must be one of: L, mL, µL/uL.")
    return u

def conc_to_M(value: float, unit_norm: str) -> float:
    return value * CONC_FACTORS[unit_norm]

def vol_to_L(value: float, unit_norm: str) -> float:
    return value * VOL_FACTORS[unit_norm]

def L_to_unit(value_L: float, unit_norm: str) -> float:
    return value_L / VOL_FACTORS[unit_norm]

def unit_lowercase_pretty(unit_norm: str) -> str:
    return {"L": "l", "ML": "ml", "UL": "ul"}[unit_norm]

def main():
    ap = argparse.ArgumentParser(description="Dilution Calculator (JSON only): compute V1 using C1·V1 = C2·V2.")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--prefix", default="dilution", help="Output JSON prefix")

    # Stock concentration
    ap.add_argument("--stock-conc", type=_float, required=True, help="Stock concentration (numeric)")
    ap.add_argument("--stock-unit", required=True, help="Stock concentration unit: M, mM, µM/uM, nM, pM")

    # Final target
    ap.add_argument("--final-conc", type=_float, required=True, help="Final concentration (numeric)")
    ap.add_argument("--final-unit", required=True, help="Final concentration unit: M, mM, µM/uM, nM, pM")
    ap.add_argument("--final-vol", type=_float, required=True, help="Final volume (numeric)")
    ap.add_argument("--final-vol-unit", required=True, help="Final volume unit: L, mL, µL/uL")

    # New: output unit for stock volume (optional)
    ap.add_argument("--stock-vol-unit-out", required=False, help="Output unit for V1: L, mL, µL/uL (default: same as final-vol-unit)")

    args = ap.parse_args()
    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)

    payload: Dict = {
        "tool": "dilution_calculator",
        "inputs": {
            "stock_concentration": args.stock_conc,
            "stock_unit": args.stock_unit,
            "final_concentration": args.final_conc,
            "final_unit": args.final_unit,
            "final_volume": args.final_vol,
            "final_volume_unit": args.final_vol_unit,
            "stock_vol_unit_out": args.stock_vol_unit_out or "(same as final volume unit)"
        },
        "formula": "C1·V1 = C2·V2 ⇒ V1 = (C2·V2)/C1",
    }

    try:
        if args.stock_conc <= 0:
            raise ValueError("Stock concentration must be > 0.")
        if args.final_conc <= 0:
            raise ValueError("Final concentration must be > 0.")
        if args.final_vol <= 0:
            raise ValueError("Final volume must be > 0.")

        uC1 = _norm_conc_unit(args.stock_unit)
        uC2 = _norm_conc_unit(args.final_unit)
        uV2 = _norm_vol_unit(args.final_vol_unit)

        # Choose output unit for V1 (default = same as final volume unit)
        if args.stock_vol_unit_out:
            uV1_out = _norm_vol_unit(args.stock_vol_unit_out)
        else:
            uV1_out = uV2

        C1_M = conc_to_M(args.stock_conc, uC1)
        C2_M = conc_to_M(args.final_conc, uC2)
        V2_L = vol_to_L(args.final_vol, uV2)

        if C1_M <= C2_M:
            raise ValueError("Stock concentration must be greater than final concentration to perform a dilution.")

        # Core calc in liters
        V1_L = (C2_M * V2_L) / C1_M

        # Express V1 in desired output unit (uV1_out)
        V1_out = L_to_unit(V1_L, uV1_out)
        unit_out_str = unit_lowercase_pretty(uV1_out)

        payload["normalized"] = {"C1_M": C1_M, "C2_M": C2_M, "V2_L": V2_L}
        payload["result"] = {
            "results": f"{V1_out} {unit_out_str}",
            "message": f"add {V1_out} {unit_out_str} of stock solution"
        }

    except Exception as e:
        payload["error"] = str(e)

    (out_dir / f"{args.prefix}.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    sys.exit(0)

if __name__ == "__main__":
    main()
