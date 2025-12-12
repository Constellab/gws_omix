#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, subprocess, sys, os
from pathlib import Path
from typing import List
os.environ.setdefault("MPLBACKEND", "Agg")

from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from collections import Counter
from PIL import Image

DPI_DEFAULT       = 300   # baked-in DPI (was a parameter)
SHOW_INDEX        = True  # always show position counters (was a parameter)
# ===================== Layout (narrow portrait look) ===============
PAGE_COLS_DEFAULT = 900   # columns per page (smaller => more pages)
WRAP_COLS         = 150   # columns shown per wrapped row (smaller => more lines)
MAX_PX_WIDTH      = 1200  # PNG width cap (pixels) after saving
LONG_THRESHOLD    = 400   # paginate if alignment length > max(LONG_THRESHOLD, PAGE_COLS_DEFAULT)

def log_init(log_path: Path):
    def _log(msg: str):
        print(msg, flush=True)
        try:
            with log_path.open("a") as L:
                L.write(msg + "\n")
        except Exception:
            pass
    return _log

def read_fasta_records(path: Path) -> List:
    try:
        return list(SeqIO.parse(str(path), "fasta"))
    except Exception:
        return []

def is_aligned(records: List) -> bool:
    if len(records) < 2:
        return False
    lengths = {len(str(r.seq)) for r in records}
    gaps = any("-" in str(r.seq) for r in records)
    return (len(lengths) == 1) or gaps

def write_fasta(records: List, out: Path):
    with out.open("wt", newline="") as w:
        SeqIO.write(records, w, "fasta")

def ensure_nonempty_fasta(fp: Path, log) -> bool:
    recs = read_fasta_records(fp)
    if len(recs) >= 2:
        return True
    log(f"[ERROR] FASTA has <2 records: {fp.name}")
    try:
        head = "".join(fp.read_text(errors="ignore").splitlines(True)[:20])
        log("[DEBUG] FASTA head:\n" + head)
    except Exception:
        pass
    return False

def run_mafft(in_fasta: Path, threads: int, log) -> str:
    # Mode précision : L-INS-i
    # mafft --localpair --maxiterate 1000 input > output
    cmd = [
       "mafft",
       "--localpair",
       "--maxiterate", "1000",
       "--thread", str(max(1, int(threads))),
       str(in_fasta)
    ]
    log("[CMD] " + " ".join(cmd))
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except FileNotFoundError:
        raise RuntimeError("MAFFT not found in environment.")
    if res.returncode != 0:
        raise RuntimeError((res.stderr or "").strip() or (res.stdout or "MAFFT failed").strip())
    out = (res.stdout or "").strip()
    if not out:
        raise RuntimeError("MAFFT returned empty alignment stdout.")
    return out

def slice_alignment(aln: MultipleSeqAlignment, start: int, end: int) -> MultipleSeqAlignment:
    return MultipleSeqAlignment([rec[start:end] for rec in aln])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_path",  required=True, help="Input FASTA")
    ap.add_argument("--out", dest="out_dir", required=True, help="Output folder")
    ap.add_argument("--prefix", default="msa", help="Output name prefix")
    ap.add_argument("--max-seqs", type=int, default=500, help="Max sequences to keep")
    ap.add_argument("--threads",  type=int, default=8, help="MAFFT threads")
    ap.add_argument("--align-if-needed", type=int, default=1, help="1: run MAFFT if not aligned")
    args = ap.parse_args()

    in_path = Path(args.in_path)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / f"{args.prefix}.log"
    log = log_init(log_path)

    # 1) read
    records = read_fasta_records(in_path)
    log(f"[INFO] Read {len(records)} sequences from {in_path.name}")
    if len(records) < 2:
        log("[ERROR] Need at least 2 sequences.")
        sys.exit(1)

    if len(records) > args.max_seqs:
        log(f"[INFO] Truncating to first {args.max_seqs} (was {len(records)}).")
        records = records[:args.max_seqs]
        write_fasta(records, out_dir / f"{args.prefix}.trimmed.fa")

    # 2) align -> *.aln.fa
    aligned_fa = out_dir / f"{args.prefix}.aln.fa"
    if is_aligned(records):
        write_fasta(records, aligned_fa)
        log("[INFO] Input appears aligned; written as-is.")
        if not ensure_nonempty_fasta(aligned_fa, log):
            sys.exit(3)
    else:
        if not args.align_if_needed:
            log("[ERROR] Not aligned and --align-if-needed=0.")
            sys.exit(1)
        try:
            aligned_text = run_mafft(in_path, args.threads, log)
        except Exception as e:
            log(f"[ERROR] MAFFT failed: {e}")
            sys.exit(2)
        # IMPORTANT : pas de newline=""
        aligned_fa.write_text(aligned_text + "\n")
        if not ensure_nonempty_fasta(aligned_fa, log):
            sys.exit(3)
        log("[INFO] MAFFT alignment completed (L-INS-i, maxiterate=1000).")

    # 3) load alignment
    try:
        aln: MultipleSeqAlignment = AlignIO.read(str(aligned_fa), "fasta")
    except Exception as e:
        log(f"[ERROR] Failed to read aligned FASTA: {e}")
        sys.exit(1)
    aln_len = aln.get_alignment_length()
    if aln_len == 0:
        log("[ERROR] Empty alignment length.")
        sys.exit(1)

    # 4) auto color scheme (DNA vs protein)
    letters = Counter("".join(str(rec.seq).upper() for rec in aln))
    total_letters = sum(v for k, v in letters.items() if k.isalpha())
    dna_like = 0.0
    if total_letters:
        dna_like = sum(letters.get(x, 0) for x in ["A","C","G","T","U","N"]) / total_letters
    scheme = "Nucleotide" if dna_like >= 0.70 else "Clustal"
    log(f"[INFO] Auto color scheme: {scheme} (DNA-like={dna_like:.2f})")

    # 5) paginate if long
    page_cols = PAGE_COLS_DEFAULT
    do_pages = aln_len > max(LONG_THRESHOLD, page_cols)
    if do_pages:
        log(f"[INFO] Rendering multiple pages. (len={aln_len}, page_cols={page_cols})")
    else:
        log(f"[INFO] Rendering 1 page(s). (len={aln_len}, page_cols={page_cols})")

    # 6) render with pyMSAviz and post-cap width
    try:
        from pymsaviz import MsaViz
    except Exception as e:
        log(f"[ERROR] pyMSAviz import failed: {e}")
        sys.exit(3)

    def render_png(aln_file: Path, out_png: Path):
        try:
            mv = MsaViz(
                str(aln_file),
                color_scheme=scheme,
                wrap_length=WRAP_COLS,
                show_count=SHOW_INDEX,
            )
            mv.savefig(str(out_png), dpi=DPI_DEFAULT)
            # post-shrink if needed (keeps width stable)
            try:
                im = Image.open(out_png)
                if im.width > MAX_PX_WIDTH:
                    new_h = int(im.height * (MAX_PX_WIDTH / im.width))
                    im = im.resize((MAX_PX_WIDTH, new_h), Image.LANCZOS)
                    im.save(out_png)
            except Exception as e:
                log(f"[WARN] PNG resize skipped: {e}")
            log(f"[INFO] Wrote {out_png.name}")
            return True
        except Exception as e:
            log(f"[ERROR] Rendering failed: {e}")
            return False

    ok = True
    if do_pages:
        starts = list(range(0, aln_len, page_cols))
        if starts and starts[-1] < aln_len and aln_len - starts[-1] < int(page_cols * 0.25):
            starts = starts[:-1]  # avoid tiny last page
        if not starts or starts[0] != 0:
            starts = [0] + starts
        pages = 0
        for s in starts:
            e = min(aln_len, s + page_cols)
            if e <= s:
                continue
            pages += 1
            sub = slice_alignment(aln, s, e)
            sub_fa = out_dir / f"{args.prefix}.page_{pages:02d}_{s+1}-{e}.fa"
            with sub_fa.open("wt", newline="") as w:
                SeqIO.write(list(sub), w, "fasta")
            png = out_dir / f"{args.prefix}.page_{pages:02d}_{s+1}-{e}.png"
            if not render_png(sub_fa, png):
                ok = False
                break
        log(f"[INFO] Total pages: {pages}")
    else:
        png = out_dir / f"{args.prefix}.png"
        if not render_png(aligned_fa, png):
            ok = False

    if not ok:
        sys.exit(3)
    log("[INFO] Done.")
    sys.exit(0)

if __name__ == "__main__":
    main()
