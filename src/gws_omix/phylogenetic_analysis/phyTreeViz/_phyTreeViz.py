#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Wrapper script around phyTreeViz CLI to generate a publication-quality
phylogenetic tree figure from a Newick (or other supported) tree file.

Intended usage (called by a Task):

  --tree           <input tree file>
  --out            <output directory>
  --prefix         <output file prefix>
  --fig-height     <figure height per leaf node>
  --fig-width      <figure width>
  --leaf-label-size <leaf label font size>
  --ignore-branch-length  (flag)
  --align-leaf-label      (flag)
  --show-branch-length    (flag)
  --show-confidence       (flag)
  --dpi            <figure dpi>

This script runs the phyTreeViz CLI roughly as:

  phytreeviz -i <tree> -o <out/prefix>.tree.png \
      --format newick \
      --fig_height <...> --fig_width <...> --leaf_label_size <...> \
      [--ignore_branch_length] [--align_leaf_label] \
      [--show_branch_length] [--show_confidence] --dpi <...>

and checks that the PNG was produced.
"""

import argparse
import subprocess
import sys
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")


def log_init(log_path: Path):
    """Create a logger that prints to stdout and appends to a log file."""
    def _log(msg: str):
        print(msg, flush=True)
        try:
            with log_path.open("a") as L:
                L.write(msg + "\n")
        except Exception:
            # Never crash because of logging
            pass
    return _log


def run_phytreeviz(
    tree_file: Path,
    out_dir: Path,
    prefix: str,
    fig_height: float,
    fig_width: float,
    leaf_label_size: int,
    ignore_branch_length: bool,
    align_leaf_label: bool,
    show_branch_length: bool,
    show_confidence: bool,
    dpi: int,
    log,
):
    """
    Run the phyTreeViz CLI with the provided options.

    Produces:
      <prefix>.tree.png              : tree figure
      <prefix>.phytreeviz_wrapper.log: wrapper log (this script)
    """
    out_png = out_dir / f"{prefix}.tree.png"

    cmd = [
        "phytreeviz",
        "-i", str(tree_file),
        "-o", str(out_png),
        "--format", "newick",  # IQ-TREE *.treefile is Newick
        "--fig_height", str(fig_height),
        "--fig_width", str(fig_width),
        "--leaf_label_size", str(leaf_label_size),
        "--dpi", str(dpi),
    ]

    if ignore_branch_length:
        cmd.append("--ignore_branch_length")
    if align_leaf_label:
        cmd.append("--align_leaf_label")
    if show_branch_length:
        cmd.append("--show_branch_length")
    if show_confidence:
        cmd.append("--show_confidence")

    log("[CMD] " + " ".join(cmd))

    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except FileNotFoundError:
        raise RuntimeError(
            "phytreeviz CLI not found in PATH. Make sure the 'phytreeviz' "
            "command is installed in this environment."
        )

    if res.returncode != 0:
        stderr = (res.stderr or "").strip()
        stdout = (res.stdout or "").strip()
        if stderr:
            log("[PHYTreeViz STDERR]\n" + stderr)
        if stdout:
            log("[PHYTreeViz STDOUT]\n" + stdout)
        raise RuntimeError(
            f"phyTreeViz CLI failed with exit code {res.returncode}. "
            f"See messages above for details."
        )

    if not out_png.is_file():
        # List files in out_dir for debugging
        files = "\n".join(sorted(p.name for p in out_dir.iterdir()))
        log(f"[ERROR] Output PNG not found: {out_png}")
        log("[ERROR] Files in out_dir:\n" + files)
        raise RuntimeError("phyTreeViz reported success but no figure was found.")


def main():
    ap = argparse.ArgumentParser(
        description="Wrapper around phyTreeViz CLI to render a phylogenetic tree figure."
    )
    ap.add_argument("--tree", dest="tree_path", required=True,
                    help="Input phylogenetic tree file (e.g. Newick from IQ-TREE).")
    ap.add_argument("--out", dest="out_dir", required=True,
                    help="Output directory for figure.")
    ap.add_argument("--prefix", default="tree",
                    help="Output prefix for the figure (default: 'tree').")

    # Figure appearance options (advanced)
    ap.add_argument(
        "--fig-height",
        type=float,
        default=0.3,
        help="Figure height per leaf node of the tree (Default: 0.3)."
    )
    ap.add_argument(
        "--fig-width",
        type=float,
        default=12.0,
        help="Figure width (Default: 12.0)."
    )
    ap.add_argument(
        "--leaf-label-size",
        type=int,
        default=8,
        help="Leaf label font size (Default: 8)."
    )
    ap.add_argument(
        "--ignore-branch-length",
        action="store_true",
        help="Ignore branch lengths when plotting (Default: use branch lengths)."
    )
    ap.add_argument(
        "--align-leaf-label",
        action="store_true",
        help="Align leaf label positions vertically (Default: OFF)."
    )
    ap.add_argument(
        "--show-branch-length",
        action="store_true",
        help="Show branch length values on the tree (Default: OFF)."
    )
    ap.add_argument(
        "--show-confidence",
        action="store_true",
        help="Show confidence / bootstrap values on the tree (Default: OFF)."
    )
    ap.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Figure DPI (Default: 300)."
    )

    args = ap.parse_args()

    tree_path = Path(args.tree_path)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    log_path = out_dir / f"{args.prefix}.phytreeviz_wrapper.log"
    log = log_init(log_path)

    log(f"[INFO] Input tree: {tree_path}")
    log(f"[INFO] Output dir: {out_dir}")
    log(f"[INFO] Prefix     : {args.prefix}")
    log(f"[INFO] fig_height : {args.fig_height}")
    log(f"[INFO] fig_width  : {args.fig_width}")
    log(f"[INFO] leaf_label_size : {args.leaf_label_size}")
    log(f"[INFO] ignore_branch_length: {args.ignore_branch_length}")
    log(f"[INFO] align_leaf_label    : {args.align_leaf_label}")
    log(f"[INFO] show_branch_length  : {args.show_branch_length}")
    log(f"[INFO] show_confidence     : {args.show_confidence}")
    log(f"[INFO] dpi                 : {args.dpi}")

    if not tree_path.is_file():
        log(f"[ERROR] Input tree file does not exist: {tree_path}")
        sys.exit(1)

    try:
        run_phytreeviz(
            tree_file=tree_path,
            out_dir=out_dir,
            prefix=args.prefix,
            fig_height=args.fig_height,
            fig_width=args.fig_width,
            leaf_label_size=args.leaf_label_size,
            ignore_branch_length=args.ignore_branch_length,
            align_leaf_label=args.align_leaf_label,
            show_branch_length=args.show_branch_length,
            show_confidence=args.show_confidence,
            dpi=args.dpi,
            log=log,
        )
    except Exception as e:
        log(f"[ERROR] phyTreeViz wrapper failed: {e}")
        sys.exit(2)

    fig_path = out_dir / f"{args.prefix}.tree.png"
    if not fig_path.is_file():
        files = "\n".join(sorted(p.name for p in out_dir.iterdir()))
        log(f"[ERROR] Expected figure not found: {fig_path}")
        log("[ERROR] Files in out_dir:\n" + files)
        sys.exit(3)

    log(f"[INFO] Figure written: {fig_path}")
    log("[INFO] phyTreeViz wrapper done.")
    sys.exit(0)


if __name__ == "__main__":
    main()
