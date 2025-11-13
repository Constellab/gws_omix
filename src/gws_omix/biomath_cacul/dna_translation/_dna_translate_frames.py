#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DNA → Protein (Forward frames 1–3) + ORFs — Text report + per-frame TSVs
========================================================================

- Reads a DNA sequence (FASTA or plain text; only the first sequence is used).
- Translates FORWARD frames 1, 2, 3 with the Standard genetic code (stops shown as '*').
- Detects ORFs in each frame: ATG → nearest in-frame stop (TAA/TAG/TGA).
- Writes:
  • <prefix>.txt  (human-readable report)  <-- NOW: ONLY RAW TRANSLATIONS, NO ORF LISTING
  • <prefix>.frame1_orfs.tsv
    <prefix>.frame2_orfs.tsv
    <prefix>.frame3_orfs.tsv   (one table per frame, first row is the header)

Dependencies:
  Biopython (Bio.Seq) — pip install biopython
"""

import argparse, sys, re, csv
from pathlib import Path
from typing import List, Tuple

try:
    from Bio.Seq import Seq
except Exception:
    print("[ERROR] Biopython is required. Install with: pip install biopython", file=sys.stderr)
    sys.exit(2)

IUPAC_DNA = set("ACGTRYSWKMBDHVNacgtryswkmbdhvn")
STANDARD_CODE = 1
STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODON = "ATG"


def _read_first_sequence(fp: Path) -> str:
    text = fp.read_text(encoding="utf-8", errors="ignore")
    lines = [ln.rstrip("\n\r") for ln in text.splitlines()]
    seq_lines: List[str] = []
    if any(ln.startswith(">") for ln in lines):
        in_seq = False
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
            if ln and not ln.startswith(">"):
                seq_lines.append(ln.strip())
    seq = re.sub(r"[\s\t\r\n]", "", "".join(seq_lines))
    if not seq:
        raise ValueError("No sequence content found.")
    return seq


def _validate_dna(seq: str):
    bad = sorted(set(ch for ch in seq if ch not in IUPAC_DNA))
    if bad:
        raise ValueError(f"Invalid DNA characters: {bad}. Allowed IUPAC DNA letters (including N).")


def _translate_frame(dna: str, start_offset: int) -> str:
    sub = dna[start_offset:]
    usable_len = (len(sub) // 3) * 3
    if usable_len <= 0:
        return ""
    sub = sub[:usable_len]
    prot = Seq(sub).translate(table=STANDARD_CODE, to_stop=False, cds=False)
    return str(prot)


def _wrap70(s: str) -> str:
    return "\n".join(s[i:i+70] for i in range(0, len(s), 70))


def _find_orfs_in_frame(dna: str, frame_offset: int) -> List[Tuple[int,int,int,int,int,str,str]]:
    """
    ORFs in a given forward frame (offset 0/1/2).
    Returns list of tuples:
      (frame, start_nt_1b, end_nt_1b, length_nt, length_aa, stop_codon, peptide)
    """
    orfs: List[Tuple[int,int,int,int,int,str,str]] = []
    n = len(dna)
    i = frame_offset
    while i + 3 <= n:
        codon = dna[i:i+3]
        if codon == START_CODON:
            j = i + 3
            while j + 3 <= n:
                c = dna[j:j+3]
                if c in STOP_CODONS:
                    start_nt_1b = i + 1
                    end_nt_1b = j + 3
                    length_nt = end_nt_1b - start_nt_1b + 1
                    cds_nt = dna[i:j]  # exclude STOP
                    pep = str(Seq(cds_nt).translate(table=STANDARD_CODE, to_stop=False, cds=False))
                    aa_len = len(pep)
                    frame = {0:1, 1:2, 2:3}[frame_offset]
                    orfs.append((frame, start_nt_1b, end_nt_1b, length_nt, aa_len, c, pep))
                    break
                j += 3
            i += 3
            continue
        i += 3

    orfs.sort(key=lambda x: (x[1], x[2]))
    return orfs


def _write_orf_tsv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        # Header on FIRST LINE:
        w.writerow(["frame", "start_nt", "end_nt", "len_nt", "len_aa", "stop", "peptide"])
        for (frame, s1, e1, ln_nt, ln_aa, stop_c, pep) in rows:
            w.writerow([frame, s1, e1, ln_nt, ln_aa, stop_c, pep])


def main():
    ap = argparse.ArgumentParser(description="Translate DNA (forward frames 1–3) + ORFs → report + TSVs")
    ap.add_argument("--in", dest="in_path", required=True, help="Input DNA file (FASTA or plain).")
    ap.add_argument("--out", required=True, help="Output directory.")
    ap.add_argument("--prefix", default="dna_translate_orf", help="Output prefix (default dna_translate_orf).")
    args = ap.parse_args()

    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)
    out_txt = out_dir / f"{args.prefix}.txt"

    # per-frame TSVs
    tsv1 = out_dir / f"{args.prefix}.frame1_orfs.tsv"
    tsv2 = out_dir / f"{args.prefix}.frame2_orfs.tsv"
    tsv3 = out_dir / f"{args.prefix}.frame3_orfs.tsv"

    try:
        dna = _read_first_sequence(Path(args.in_path))
        _validate_dna(dna)
        dna_up = dna.upper()

        # Raw frame translations
        frames = []
        for (frame, off) in [(1,0),(2,1),(3,2)]:
            prot = _translate_frame(dna_up, off)
            frames.append((frame, prot, off))

        # ORFs (still computed and written to TSVs)
        orfs_f1 = _find_orfs_in_frame(dna_up, 0)
        orfs_f2 = _find_orfs_in_frame(dna_up, 1)
        orfs_f3 = _find_orfs_in_frame(dna_up, 2)

        # Write per-frame TSVs (header on first line)
        _write_orf_tsv(tsv1, orfs_f1)
        _write_orf_tsv(tsv2, orfs_f2)
        _write_orf_tsv(tsv3, orfs_f3)

        # Report (NOW ONLY RAW TRANSLATIONS, NO ORF LISTING)
        with out_txt.open("w", encoding="utf-8") as w:
            w.write("DNA → Protein (Forward frames 1–3)\n")
            w.write("Standard genetic code. Stops shown as '*'.\n")
            w.write(f"Input length: {len(dna_up)} nt\n\n")

            w.write("=== Raw frame translations (forward) ===\n\n")
            for frame, prot, off in frames:
                w.write(f">protein sequence, frame {frame} (strand +, offset {off})\n")
                w.write((_wrap70(prot) or "(empty)") + "\n\n")

        print(f"[OK] Wrote report: {out_txt}")
        print(f"[OK] Wrote TSVs: {tsv1.name}, {tsv2.name}, {tsv3.name}")

    except Exception as e:
        with out_txt.open("w", encoding="utf-8") as w:
            w.write("ERROR\n")
            w.write(str(e) + "\n")
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(0)


if __name__ == "__main__":
    main()
