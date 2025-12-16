#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import subprocess
import sys
import os
from pathlib import Path


def log_init(log_path: Path):
    """Create a simple logging function that writes both to stdout and a log file."""
    def _log(msg: str):
        print(msg, flush=True)
        try:
            with log_path.open("a") as L:
                L.write(msg + "\n")
        except Exception:
            pass
    return _log


def run_iqtree(
    in_fasta: Path,
    out_dir: Path,
    prefix: str,
    model: str,
    bootstrap: int,
    threads: int,
    seed: int,
    log,
):
    """
    Run IQ-TREE on an aligned FASTA file.

    Produces at minimum:
      <prefix>.treefile   (Newick tree)
      <prefix>.log        (IQ-TREE log)
      <prefix>.iqtree     (run summary)

    All outputs are written with -pre <out_dir>/<prefix>.
    """

    pre = out_dir / prefix

    cmd = [
        "iqtree",                      # force use of 'iqtree' binary
        "-s", str(in_fasta),
        "-m", model,
        "-nt", str(max(1, int(threads))),
        "-pre", str(pre),
    ]

    # Ultrafast bootstrap
    if bootstrap and bootstrap > 0:
        cmd += ["-bb", str(int(bootstrap))]

    # Optional seed (0 = IQ-TREE default)
    if seed and seed > 0:
        cmd += ["-seed", str(int(seed))]

    log("[CMD] " + " ".join(cmd))

    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except FileNotFoundError:
        raise RuntimeError(
            "iqtree not found in PATH. Make sure the 'iqtree' binary is "
            "installed and available in this environment."
        )

    if res.returncode != 0:
        stderr = (res.stderr or "").strip()
        stdout = (res.stdout or "").strip()
        if stderr:
            log("[IQTREE STDERR]\n" + stderr)
        if stdout:
            log("[IQTREE STDOUT]\n" + stdout)
        raise RuntimeError(
            f"IQ-TREE (iqtree) failed with exit code {res.returncode}. "
            f"See IQ-TREE log file and messages above."
        )


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Wrapper around IQ-TREE to infer a phylogenetic tree "
            "from an aligned FASTA file."
        )
    )
    ap.add_argument("--in", dest="in_path", required=True,
                    help="Input aligned FASTA")
    ap.add_argument("--out", dest="out_dir", required=True,
                    help="Output directory")
    ap.add_argument("--prefix", default="tree",
                    help="Output prefix for IQ-TREE files")
    ap.add_argument(
        "--model",
        default="MFP",
        help="Substitution model (e.g. MFP, GTR+G, LG+G). "
             "MFP = ModelFinder Plus (automatic model selection)."
    )
    ap.add_argument(
        "--bootstrap",
        type=int,
        default=1000,
        help="Ultrafast bootstrap replicates (-bb). 0 = no bootstrap."
    )
    ap.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of CPU threads (-nt)."
    )
    ap.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed for IQ-TREE (0 = IQ-TREE default)."
    )
    args = ap.parse_args()

    in_path = Path(args.in_path)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Wrapper log (in addition to IQ-TREE's own <prefix>.log)
    # NOTE: name kept as '.iqtree2_wrapper.log' to match the Task implementation.
    log_path = out_dir / f"{args.prefix}.iqtree2_wrapper.log"
    log = log_init(log_path)

    log(f"[INFO] Input FASTA: {in_path}")
    log(f"[INFO] Output dir : {out_dir}")
    log(f"[INFO] Prefix      : {args.prefix}")
    log(f"[INFO] Model       : {args.model}")
    log(f"[INFO] Bootstrap   : {args.bootstrap}")
    log(f"[INFO] Threads     : {args.threads}")
    log(f"[INFO] Seed        : {args.seed}")

    if not in_path.is_file():
        log(f"[ERROR] Input FASTA does not exist: {in_path}")
        sys.exit(1)

    try:
        run_iqtree(
            in_fasta=in_path,
            out_dir=out_dir,
            prefix=args.prefix,
            model=args.model,
            bootstrap=args.bootstrap,
            threads=args.threads,
            seed=args.seed,
            log=log,
        )
    except Exception as e:
        log(f"[ERROR] IQ-TREE run failed: {e}")
        sys.exit(2)

    # IQ-TREE should have produced <prefix>.treefile
    treefile = out_dir / f"{args.prefix}.treefile"
    if not treefile.is_file():
        # List files in out_dir for easier debugging
        files = "\n".join(sorted(p.name for p in out_dir.iterdir()))
        log(f"[ERROR] Treefile not found: {treefile}")
        log("[ERROR] Files in out_dir:\n" + files)
        sys.exit(3)

    log(f"[INFO] Treefile: {treefile}")
    log("[INFO] IQ-TREE wrapper done.")
    sys.exit(0)


if __name__ == "__main__":
    main()
