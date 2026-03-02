#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
_gaf_to_gmt.py — Convert GO GAF to GMT gene sets.

This script always generates:
A) GO_ALL (4 GMT):
  - <prefix>_GO_ALL_withIEA_id.gmt
  - <prefix>_GO_ALL_withIEA_symbol.gmt
  - <prefix>_GO_ALL_noIEA_id.gmt
  - <prefix>_GO_ALL_noIEA_symbol.gmt

B) GO split by aspect BP/MF/CC (12 GMT):
  - <prefix>_GO_BP_withIEA_id.gmt
  - <prefix>_GO_BP_withIEA_symbol.gmt
  - <prefix>_GO_MF_withIEA_id.gmt
  - <prefix>_GO_MF_withIEA_symbol.gmt
  - <prefix>_GO_CC_withIEA_id.gmt
  - <prefix>_GO_CC_withIEA_symbol.gmt
  - <prefix>_GO_BP_noIEA_id.gmt
  - <prefix>_GO_BP_noIEA_symbol.gmt
  - <prefix>_GO_MF_noIEA_id.gmt
  - <prefix>_GO_MF_noIEA_symbol.gmt
  - <prefix>_GO_CC_noIEA_id.gmt
  - <prefix>_GO_CC_noIEA_symbol.gmt

Design choices:
  - Always produce both withIEA and noIEA.
  - NOT-qualified annotations are always excluded.
  - No term-size filtering here (filter sizes at GSEA step).
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, Set, Tuple

GO_RE = re.compile(r"^GO:\d{7}$")
ASPECTS = {"P": "BP", "F": "MF", "C": "CC"}  # GAF Aspect col: P/F/C


def eprint(*args) -> None:
    print(*args, file=sys.stderr)


def open_maybe_gz(path: Path):
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, "r", encoding="utf-8", errors="replace")


def _new_maps() -> Tuple[
    DefaultDict[str, Set[str]],
    DefaultDict[str, Set[str]],
    DefaultDict[str, Set[str]],
    DefaultDict[str, Set[str]],
]:
    with_id: DefaultDict[str, Set[str]] = defaultdict(set)
    with_sym: DefaultDict[str, Set[str]] = defaultdict(set)
    no_id: DefaultDict[str, Set[str]] = defaultdict(set)
    no_sym: DefaultDict[str, Set[str]] = defaultdict(set)
    return with_id, with_sym, no_id, no_sym


def parse_gaf_maps(gaf_path: Path) -> Tuple[
    # GO_ALL
    DefaultDict[str, Set[str]],
    DefaultDict[str, Set[str]],
    DefaultDict[str, Set[str]],
    DefaultDict[str, Set[str]],
    # By aspect: BP/MF/CC -> maps
    Dict[str, DefaultDict[str, Set[str]]],
    Dict[str, DefaultDict[str, Set[str]]],
    Dict[str, DefaultDict[str, Set[str]]],
    Dict[str, DefaultDict[str, Set[str]]],
]:
    """
    Returns:
      GO_ALL: withIEA_id, withIEA_symbol, noIEA_id, noIEA_symbol
      GO_SPLIT: dict(aspect in {BP,MF,CC} -> same 4 maps)
    """
    all_with_id, all_with_sym, all_no_id, all_no_sym = _new_maps()

    split_with_id: Dict[str, DefaultDict[str, Set[str]]] = {k: defaultdict(set) for k in ASPECTS.values()}
    split_with_sym: Dict[str, DefaultDict[str, Set[str]]] = {k: defaultdict(set) for k in ASPECTS.values()}
    split_no_id:   Dict[str, DefaultDict[str, Set[str]]] = {k: defaultdict(set) for k in ASPECTS.values()}
    split_no_sym:  Dict[str, DefaultDict[str, Set[str]]] = {k: defaultdict(set) for k in ASPECTS.values()}

    with open_maybe_gz(gaf_path) as fh:
        for line in fh:
            if not line or line.startswith("!"):
                continue
            parts = line.rstrip("\n").split("\t")
            # We need at least up to Aspect column (index 8)
            if len(parts) < 9:
                continue

            db_object_id = parts[1].strip()
            db_object_symbol = parts[2].strip()
            qualifier = parts[3].strip()
            go_id = parts[4].strip()
            evidence = parts[6].strip()
            aspect_raw = parts[8].strip()  # P/F/C

            if not GO_RE.match(go_id):
                continue

            # Always exclude NOT-qualified annotations
            if "NOT" in qualifier.split("|"):
                continue

            aspect = ASPECTS.get(aspect_raw)
            if aspect is None:
                # Unknown aspect -> ignore (keeps output clean)
                continue

            # GO_ALL withIEA
            if db_object_id:
                all_with_id[go_id].add(db_object_id)
                split_with_id[aspect][go_id].add(db_object_id)
            if db_object_symbol:
                all_with_sym[go_id].add(db_object_symbol)
                split_with_sym[aspect][go_id].add(db_object_symbol)

            # noIEA branch
            if evidence != "IEA":
                if db_object_id:
                    all_no_id[go_id].add(db_object_id)
                    split_no_id[aspect][go_id].add(db_object_id)
                if db_object_symbol:
                    all_no_sym[go_id].add(db_object_symbol)
                    split_no_sym[aspect][go_id].add(db_object_symbol)

    return (
        all_with_id, all_with_sym, all_no_id, all_no_sym,
        split_with_id, split_with_sym, split_no_id, split_no_sym
    )


def write_gmt(path: Path, term_to_genes: DefaultDict[str, Set[str]]) -> int:
    """
    GMT format:
      TERM <tab> na <tab> gene1 <tab> gene2 ...
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(path, "w", encoding="utf-8") as out:
        for term in sorted(term_to_genes.keys()):
            genes = sorted(g for g in term_to_genes[term] if g)
            if not genes:
                continue
            out.write(term + "\t" + "na" + "\t" + "\t".join(genes) + "\n")
            n += 1
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description="Convert GO GAF into GMT (GO_ALL + BP/MF/CC) for withIEA/noIEA × id/symbol.")
    ap.add_argument("--gaf", required=True, help="Input GAF file path (.gaf or .gaf.gz).")
    ap.add_argument("--outdir", required=True, help="Output directory for GMT files.")
    ap.add_argument("--prefix", default="organism", help="Prefix for output files (default: organism).")
    args = ap.parse_args()

    gaf_path = Path(args.gaf)
    if not gaf_path.exists():
        eprint(f"ERROR: GAF file not found: {gaf_path}")
        return 2

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    (
        all_with_id, all_with_sym, all_no_id, all_no_sym,
        split_with_id, split_with_sym, split_no_id, split_no_sym
    ) = parse_gaf_maps(gaf_path)

    if not all_with_id and not all_with_sym:
        eprint("ERROR: No annotations parsed from GAF.")
        return 3

    # --- GO_ALL (4) ---
    out_with_id = outdir / f"{args.prefix}_GO_ALL_withIEA_id.gmt"
    out_with_sym = outdir / f"{args.prefix}_GO_ALL_withIEA_symbol.gmt"
    out_no_id = outdir / f"{args.prefix}_GO_ALL_noIEA_id.gmt"
    out_no_sym = outdir / f"{args.prefix}_GO_ALL_noIEA_symbol.gmt"

    n1 = write_gmt(out_with_id, all_with_id)
    n2 = write_gmt(out_with_sym, all_with_sym)
    n3 = write_gmt(out_no_id, all_no_id)
    n4 = write_gmt(out_no_sym, all_no_sym)

    eprint(f"[INFO] wrote {n1} terms -> {out_with_id.name}")
    eprint(f"[INFO] wrote {n2} terms -> {out_with_sym.name}")
    eprint(f"[INFO] wrote {n3} terms -> {out_no_id.name}")
    eprint(f"[INFO] wrote {n4} terms -> {out_no_sym.name}")

    # --- GO split by aspect (12) ---
    for aspect in ("BP", "MF", "CC"):
        fp1 = outdir / f"{args.prefix}_GO_{aspect}_withIEA_id.gmt"
        fp2 = outdir / f"{args.prefix}_GO_{aspect}_withIEA_symbol.gmt"
        fp3 = outdir / f"{args.prefix}_GO_{aspect}_noIEA_id.gmt"
        fp4 = outdir / f"{args.prefix}_GO_{aspect}_noIEA_symbol.gmt"

        m1 = write_gmt(fp1, split_with_id[aspect])
        m2 = write_gmt(fp2, split_with_sym[aspect])
        m3 = write_gmt(fp3, split_no_id[aspect])
        m4 = write_gmt(fp4, split_no_sym[aspect])

        eprint(f"[INFO] wrote {m1} terms -> {fp1.name}")
        eprint(f"[INFO] wrote {m2} terms -> {fp2.name}")
        eprint(f"[INFO] wrote {m3} terms -> {fp3.name}")
        eprint(f"[INFO] wrote {m4} terms -> {fp4.name}")

    eprint("✔ Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())