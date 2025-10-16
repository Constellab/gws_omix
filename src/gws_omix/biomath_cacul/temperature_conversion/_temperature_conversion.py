#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, json, sys
from pathlib import Path
from typing import Dict, Any


def convert_from_celsius(c: float) -> Dict[str, float]:
    f = (c * 9 / 5) + 32
    k = c + 273.16
    return {"C": c, "F": f, "K": k}


def convert_from_fahrenheit(f: float) -> Dict[str, float]:
    c = (f - 32) * 5 / 9
    k = c + 273.16
    return {"C": c, "F": f, "K": k}


def convert_from_kelvin(k: float) -> Dict[str, float]:
    c = k - 273.16
    f = (c * 9 / 5) + 32
    return {"C": c, "F": f, "K": k}


def main():
    ap = argparse.ArgumentParser(description="Temperature Conversion (°C, °F, °K)")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--prefix", default="temperature_conversion", help="Output JSON prefix")
    ap.add_argument("--value", type=float, required=True, help="Temperature value to convert")
    ap.add_argument("--unit", required=True, choices=["C", "F", "K"], help="Input temperature unit (C, F, or K)")
    args = ap.parse_args()

    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)

    payload: Dict[str, Any] = {
        "tool": "temperature_conversion",
        "input": {"value": args.value, "unit": args.unit},
        "formula": {
            "C": "°C = (°F − 32) × 5/9",
            "F": "°F = (°C × 9/5) + 32",
            "K": "°K = °C + 273.16"
        }
    }

    try:
        if args.unit == "C":
            result = convert_from_celsius(args.value)
        elif args.unit == "F":
            result = convert_from_fahrenheit(args.value)
        elif args.unit == "K":
            result = convert_from_kelvin(args.value)
        else:
            raise ValueError("Unit must be one of: C, F, K.")

        payload["result"] = {
            "°C": round(result["C"], 2),
            "°F": round(result["F"], 2),
            "°K": round(result["K"], 2),
        }

    except Exception as e:
        payload["error"] = str(e)

    (out_dir / f"{args.prefix}.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    sys.exit(0)


if __name__ == "__main__":
    main()
