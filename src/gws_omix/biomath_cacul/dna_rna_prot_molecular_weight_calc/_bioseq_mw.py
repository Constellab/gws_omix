#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse, json, sys, re
from pathlib import Path
from typing import Dict, Any, List

# ---- DNA masses (vendor-precision; matches vendor calculators) ----
DNA_MASS = {"A": 313.209, "T": 304.197, "C": 289.184, "G": 329.205}

# ---- RNA masses (vendor-precision) ----
RNA_MASS = {"A": 329.21, "U": 306.17, "C": 305.18, "G": 345.21}

# ---- dsDNA / dsRNA base-pair averages (Da per bp) tunées ----
DS_DNA_BP_AVG = 618.004812
DS_RNA_BP_AVG = 644.574228

# ---- Perte d'eau (Da) lors de la ligation circulaire ----
H2O = 18.01528

# ---- Deltas d'extrémités (Da) pour coller aux sorties du site ----
# DNA (par brin)
SS_DNA_END_HYDROXYL     = -67.497
SS_DNA_END_PHOSPHATE    = SS_DNA_END_HYDROXYL + 79.98
SS_DNA_END_TRIPHOSPHATE = SS_DNA_END_HYDROXYL + 239.94

# dsDNA (delta total molécule)
DS_DNA_END = {
    "hydroxyl":   -123.92,
    "phosphate":   +36.04,
    "triphosphate":+355.96,
}

# RNA (par brin)
SS_RNA_END_HYDROXYL     = -72.10
SS_RNA_END_PHOSPHATE    = SS_RNA_END_HYDROXYL + 79.98
SS_RNA_END_TRIPHOSPHATE = SS_RNA_END_HYDROXYL + 159.96

TYPE_CHOICES = {"DNA", "RNA", "PROTEIN"}
DNA_ALLOWED = set("ATGCatgc")
RNA_ALLOWED = set("AUGCaugc")
AA_ALLOWED_20 = set("ACDEFGHIKLMNPQRSTVWY")  # 20 AAs pour ProtParam

def _clean(s: str) -> str:
    return re.sub(r"[\s\r\n\t]", "", s or "")

def _read_first_sequence(fp: Path) -> str:
    text = fp.read_text(encoding="utf-8", errors="ignore")
    lines = [ln.rstrip("\n\r") for ln in text.splitlines()]
    seq_lines: List[str] = []
    in_seq = False
    saw_fasta = any(ln.startswith(">") for ln in lines)
    if saw_fasta:
        for ln in lines:
            if ln.startswith(">"):
                if in_seq:
                    break
                in_seq = True
                continue
            if in_seq:
                seq_lines.append(ln.strip())
    else:
        for ln in lines:
            if ln.startswith(">"):
                continue
            if ln.strip():
                seq_lines.append(ln.strip())
    seq = _clean("".join(seq_lines))
    if not seq:
        raise ValueError("No sequence content found in input file.")
    return seq

def _validate_na(seq: str, na_type: str):
    s = seq.upper()
    if not s:
        raise ValueError("Empty nucleic acid sequence.")
    if na_type == "DNA":
        bad = sorted(set(ch for ch in s if ch not in DNA_ALLOWED))
        if bad:
            raise ValueError(f"Invalid DNA characters: {bad}. Allowed: A/T/G/C only.")
    else:
        bad = sorted(set(ch for ch in s if ch not in RNA_ALLOWED))
        if bad:
            raise ValueError(f"Invalid RNA characters: {bad}. Allowed: A/U/G/C only.")

def _na_counts(seq: str, na_type: str) -> Dict[str, int]:
    up = seq.upper()
    keys = ("A", "T", "C", "G") if na_type == "DNA" else ("A", "U", "C", "G")
    return {k: up.count(k) for k in keys}

def _na_sum_mass(counts: Dict[str, int], table: Dict[str, float]) -> float:
    return sum(counts[b] * table[b] for b in counts.keys())

def _per_strand_delta(na_type: str, fivep: str) -> float:
    f = (fivep or "").lower()
    if na_type == "DNA":
        if f == "hydroxyl":     return SS_DNA_END_HYDROXYL
        if f == "phosphate":    return SS_DNA_END_PHOSPHATE
        if f == "triphosphate": return SS_DNA_END_TRIPHOSPHATE
        raise ValueError("five_prime must be one of: hydroxyl, phosphate, triphosphate")
    else:  # RNA
        if f == "hydroxyl":     return SS_RNA_END_HYDROXYL
        if f == "phosphate":    return SS_RNA_END_PHOSPHATE
        if f == "triphosphate": return SS_RNA_END_TRIPHOSPHATE
        raise ValueError("five_prime must be one of: hydroxyl, phosphate, triphosphate")

# ------------------------- PROTEIN (ProtParam) -------------------------

def _validate_protein(seq: str):
    s = seq.upper()
    if not s:
        raise ValueError("Empty protein sequence.")
    bad = sorted(set(ch for ch in s if ch not in AA_ALLOWED_20))
    if bad:
        raise ValueError(f"Unsupported amino acids: {bad}. Allowed: {sorted(AA_ALLOWED_20)}")

def _protein_analysis(seq: str) -> Dict[str, Any]:
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
    except Exception as e:
        raise RuntimeError(
            f"Biopython not available ({e}). Install with: conda install -c conda-forge biopython"
        )

    up = seq.upper()
    pa = ProteinAnalysis(up)

    # Comptes & pourcentages (map complet sur 20 AA ; 0 pour AA absents)
    counts = {aa: up.count(aa) for aa in sorted(AA_ALLOWED_20)}
    perc   = pa.get_amino_acids_percent()
    perc_full = {aa: round(float(perc.get(aa, 0.0)) * 100, 2) for aa in sorted(AA_ALLOWED_20)}

    # Mass + autres indicateurs
    mw = float(pa.molecular_weight())  # average mass, H2O incluse
    arom = float(pa.aromaticity())
    instab = float(pa.instability_index())
    pI = float(pa.isoelectric_point())
    helix, turn, sheet = pa.secondary_structure_fraction()
    eps_red, eps_ox = pa.molar_extinction_coefficient()  # (reduced, oxidized)

    return {
        "counts": counts,
        "percent": perc_full,
        "molecular_weight_Da_raw": mw,
        "molecular_weight_Da": round(mw, 2),
        "aromaticity": round(arom, 4),
        "instability_index": round(instab, 2),
        "isoelectric_point": round(pI, 2),
        "secondary_structure_fraction": {
            "helix": round(float(helix), 4),
            "turn":  round(float(turn), 4),
            "sheet": round(float(sheet), 4),
        },
        "molar_extinction_coefficient": {
            "reduced": int(eps_red),
            "oxidized": int(eps_ox),
        }
    }

# ------------------------------ MAIN -----------------------------------

def main():
    ap = argparse.ArgumentParser(description="BioSeq MW (DNA/RNA/PROTEIN) from sequence file (JSON only).")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--prefix", default="bioseq_mw", help="Output JSON prefix")

    ap.add_argument("--type", required=True, choices=sorted(TYPE_CHOICES), help="DNA | RNA | PROTEIN")
    ap.add_argument("--in", dest="in_path", required=True, help="Path to sequence file (FASTA or plain sequence)")

    # NA-only (doit être renseigné pour DNA/RNA ; ignoré pour PROTEIN)
    ap.add_argument("--strand", default="", help="NA only: '' | ss | ds (must be non-empty for DNA/RNA)")
    ap.add_argument("--topology", default="", help="NA only: '' | linear | circular (must be non-empty for DNA/RNA)")
    ap.add_argument("--five-prime", default="", help="NA only: '' | hydroxyl | phosphate | triphosphate (must be non-empty for DNA/RNA)")

    args = ap.parse_args()
    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)
    out_json = out_dir / f"{args.prefix}.json"

    try:
        seq = _read_first_sequence(Path(args.in_path))

        # -------------------- DNA / RNA --------------------
        if args.type in ("DNA", "RNA"):
            strand   = (args.strand or "").strip().lower()
            topology = (args.topology or "").strip().lower()
            fivep_in = (args.five_prime or "").strip().lower()

            missing = []
            if strand not in ("ss", "ds"):
                missing.append("strand (ss|ds)")
            if topology not in ("linear", "circular"):
                missing.append("topology (linear|circular)")
            if fivep_in not in ("hydroxyl", "phosphate", "triphosphate"):
                missing.append("five_prime (hydroxyl|phosphate|triphosphate)")

            if missing:
                payload = {
                    "tool": "bioseq_mw",
                    "input": {
                        "type": args.type,
                        "sequence": seq,
                        "length": len(seq),
                        "strand": args.strand,
                        "topology": args.topology,
                        "five_prime": args.five_prime
                    },
                    "error": "Missing/invalid NA-only parameters: " + ", ".join(missing),
                    "hint": "For DNA/RNA, set strand=ss|ds, topology=linear|circular, five_prime=hydroxyl|phosphate|triphosphate."
                }
                out_json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
                print(json.dumps(payload, indent=2, ensure_ascii=False))
                sys.exit(0)

            # Validation
            _validate_na(seq, args.type)
            counts = _na_counts(seq, args.type)

            # Circular => force 5' phosphate pour ligation
            fivep_effective = "phosphate" if topology == "circular" else fivep_in

            if args.type == "DNA":
                if strand == "ss":
                    one_strand_sum = _na_sum_mass(counts, DNA_MASS)
                    per_strand     = _per_strand_delta("DNA", fivep_effective)
                    mw = one_strand_sum + per_strand
                    debug = {"mode": "ss", "one_strand_sum": one_strand_sum, "per_strand_delta": per_strand}
                    if topology == "circular":  # fermeture d’un brin
                        mw -= H2O
                        debug["circular_ligation_delta"] = -H2O
                else:
                    n_bp = len(seq)
                    ds_end_delta = DS_DNA_END[fivep_effective]
                    mw = n_bp * DS_DNA_BP_AVG + ds_end_delta
                    debug = {"mode": "ds_avg", "n_bp": n_bp, "ds_basepair_avg": DS_DNA_BP_AVG, "ds_end_delta": ds_end_delta}
                    if topology == "circular":  # fermeture de 2 brins
                        mw -= 2.0 * H2O
                        debug["circular_ligation_delta"] = -2.0 * H2O

            else:  # RNA
                one_strand_sum = _na_sum_mass(counts, RNA_MASS)
                per_strand     = _per_strand_delta("RNA", fivep_effective)
                if strand == "ss":
                    mw = one_strand_sum + per_strand
                    debug = {"mode": "ss", "one_strand_sum": one_strand_sum, "per_strand_delta": per_strand}
                    if topology == "circular":
                        mw -= H2O
                        debug["circular_ligation_delta"] = -H2O
                else:
                    n_bp = len(seq)
                    mw = n_bp * DS_RNA_BP_AVG
                    debug = {"mode": "ds_avg_rna", "n_bp": n_bp, "ds_basepair_avg_rna": DS_RNA_BP_AVG}
                    if topology == "circular":
                        mw -= 2.0 * H2O
                        debug["circular_ligation_delta"] = -2.0 * H2O

            mw_rounded = round(mw, 2)

            payload: Dict[str, Any] = {
                "tool": "bioseq_mw",
                "input": {
                    "type": args.type,
                    "sequence": seq,
                    "length": len(seq),
                    "strand": strand,
                    "topology": topology,
                    "five_prime": fivep_in,
                    "five_prime_effective": fivep_effective if topology == "circular" else fivep_in
                },
                "counts": counts,
                "masses_used": {"DNA": DNA_MASS if args.type == "DNA" else None,
                                "RNA": RNA_MASS if args.type == "RNA" else None},
                "formula": {
                    "DNA_ss": "Σ(A*313.209 + T*304.197 + C*289.184 + G*329.205) + tuned end delta",
                    "DNA_ds": f"MW = N_bp * {DS_DNA_BP_AVG} + tuned ds end delta",
                    "RNA_ss": "Σ(A*329.21 + U*306.17 + C*305.18 + G*345.21) + tuned end delta",
                    "RNA_ds": f"MW = N_bp * {DS_RNA_BP_AVG} (no end delta)",
                    "circular_ligation_delta": {"ss": -18.01528, "ds": -36.03056, "applies_if": "topology=circular (5′ phosphate forced)"}
                },
                "debug": debug,
                "result": {
                    "molecular_weight_Da": mw_rounded,
                    "message": f"{args.type} {strand}/{topology}, 5'={fivep_effective} → MW = {mw_rounded:.2f} Da"
                }
            }

        # -------------------- PROTEIN (ProtParam) --------------------
        else:
            _validate_protein(seq)
            info = _protein_analysis(seq)
            payload = {
                "tool": "bioseq_mw",
                "input": {
                    "type": "PROTEIN",
                    "sequence": seq,
                    "length": len(seq),
                    # on garde ces champs vides/ignorés pour PROTEIN (compat)
                    "strand": "",
                    "topology": "",
                    "five_prime": ""
                },
                "counts": info["counts"],
                "percent": info["percent"],  # proportion par AA (0–1)
                "formula": "Biopython ProtParam (average mass with terminal H2O)",
                "analysis": {
                    "aromaticity": info["aromaticity"],
                    "instability_index": info["instability_index"],
                    "isoelectric_point": info["isoelectric_point"],
                    "secondary_structure_fraction": info["secondary_structure_fraction"],
                    "molar_extinction_coefficient": info["molar_extinction_coefficient"]
                },
                "result": {
                    "molecular_weight_Da": info["molecular_weight_Da"],
                    "message": f"Protein MW = {info['molecular_weight_Da']:.2f} Da"
                },
                "note": "For PROTEIN, strand/topology/five_prime are ignored. Set them to empty ('')."
            }

    except Exception as e:
        payload = {"tool": "bioseq_mw", "error": str(e)}

    out_json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    sys.exit(0)

if __name__ == "__main__":
    main()
