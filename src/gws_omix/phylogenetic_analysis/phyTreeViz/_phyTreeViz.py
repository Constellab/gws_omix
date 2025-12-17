#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import argparse
import sys
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")


def log_init():
    """Simple logger: stdout only, no log file on disk."""
    def _log(msg: str):
        print(msg, flush=True)
    return _log


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Render a phylogenetic tree using phyTreeViz from a Newick tree file."
        )
    )
    ap.add_argument("--tree", dest="tree_path", required=True,
                    help="Input tree file (Newick, etc.)")
    ap.add_argument("--out", dest="out_dir", required=True,
                    help="Output directory")
    ap.add_argument("--prefix", default="tree",
                    help="Output prefix (default: 'tree')")

    # Advanced figure options (passed from the Task)
    ap.add_argument(
        "--fig-height",
        dest="fig_height",
        type=float,
        default=0.3,
        help="Figure height per leaf node (Default: 0.3)."
    )
    ap.add_argument(
        "--fig-width",
        dest="fig_width",
        type=float,
        default=12.0,
        help="Figure width (Default: 12.0)."
    )
    ap.add_argument(
        "--leaf-label-size",
        dest="leaf_label_size",
        type=float,
        default=8.0,
        help="Leaf label font size (Default: 8)."
    )
    ap.add_argument(
        "--align-leaf-label",
        dest="align_leaf_label",
        action="store_true",
        help="Align leaf label positions."
    )
    ap.add_argument(
        "--show-branch-length",
        dest="show_branch_length",
        action="store_true",
        help="Show branch length values on branches."
    )
    ap.add_argument(
        "--show-confidence",
        dest="show_confidence",
        action="store_true",
        help="Show confidence / bootstrap values on branches."
    )
    ap.add_argument(
        "--dpi",
        dest="dpi",
        type=int,
        default=300,
        help="Figure DPI (Default: 300)."
    )
    args = ap.parse_args()

    tree_path = Path(args.tree_path)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    log = log_init()

    log(f"[INFO] Input tree: {tree_path}")
    log(f"[INFO] Output dir: {out_dir}")
    log(f"[INFO] Prefix     : {args.prefix}")
    log(f"[INFO] fig_height : {args.fig_height}")
    log(f"[INFO] fig_width  : {args.fig_width}")
    log(f"[INFO] leaf_label_size : {args.leaf_label_size}")
    log(f"[INFO] align_leaf_label    : {args.align_leaf_label}")
    log(f"[INFO] show_branch_length  : {args.show_branch_length}")
    log(f"[INFO] show_confidence     : {args.show_confidence}")
    log(f"[INFO] dpi                 : {args.dpi}")

    if not tree_path.is_file():
        log(f"[ERROR] Input tree file does not exist: {tree_path}")
        sys.exit(1)

    try:
        from phytreeviz import TreeViz
    except ImportError as e:
        log("[ERROR] Could not import 'phytreeviz'. "
            "Make sure the package is installed in this environment.")
        log(f"[ERROR] ImportError: {e}")
        sys.exit(2)

    out_png = out_dir / f"{args.prefix}.tree.png"

    try:
        # Branch lengths are always used to set the horizontal scale.
        tv = TreeViz(
            tree_data=str(tree_path),
            format="newick",
            height=args.fig_height,
            width=args.fig_width,
            align_leaf_label=args.align_leaf_label,
            leaf_label_size=args.leaf_label_size,
        )

        # Branch-length and confidence labels (text)
        if args.show_branch_length:
            tv.show_branch_length(size=max(6, int(args.leaf_label_size * 0.7)))
        if args.show_confidence:
            tv.show_confidence(size=max(6, int(args.leaf_label_size * 0.7)))

        # Always draw a scale bar
        tv.show_scale_bar()

        tv.savefig(str(out_png), dpi=args.dpi)
        log(f"[INFO] Wrote tree image: {out_png}")

    except Exception as e:
        log(f"[ERROR] phyTreeViz rendering failed: {e}")
        sys.exit(3)

    if not out_png.is_file():
        log(f"[ERROR] Output PNG not found: {out_png}")
        sys.exit(4)

    log("[INFO] phyTreeViz wrapper done.")
    sys.exit(0)


if __name__ == "__main__":
    main()
