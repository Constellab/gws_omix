#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank


def setup_logger():
    logging.basicConfig(
        level=logging.INFO,
        format="[%(levelname)s] %(message)s",
    )


def parse_args():
    parser = argparse.ArgumentParser(description="Genome visualization using pyGenomeViz")
    parser.add_argument("--in", dest="input_dir", required=True)
    parser.add_argument("--out", dest="output_dir", required=True)
    parser.add_argument("--prefix", required=True)
    return parser.parse_args()


def main():
    setup_logger()
    args = parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    prefix = args.prefix

    output_dir.mkdir(parents=True, exist_ok=True)
    out_html = output_dir / f"{prefix}.html"

    logging.info(f"Input directory: {input_dir}")
    logging.info(f"Output directory: {output_dir}")
    logging.info(f"Prefix: {prefix}")

    gbk_files = sorted(input_dir.glob("*.gb"))
    if not gbk_files:
        raise RuntimeError("No GenBank files found")

    logging.info(f"Found {len(gbk_files)} GenBank files")

    gv = GenomeViz()
    gv.set_scale_bar(ymargin=0.5)

    track_count = 0

    for gbk_file in gbk_files:
        try:
            gbk = Genbank(gbk_file)
            logging.info(
                f"Loaded {gbk_file.name}, genome length: {gbk.genome_length}"
            )

            # ✅ SEGMENTS = genome_length (OBLIGATOIRE)
            track = gv.add_feature_track(
                gbk.name or gbk_file.stem,
                gbk.genome_length,
            )

            track.add_sublabel()

            features = gbk.extract_features(feature_type="CDS")

            if not features:
                logging.warning(f"No CDS features in {gbk_file.name}")
                continue

            track.add_features(features)
            track_count += 1

        except Exception as e:
            logging.warning(f"Failed to process {gbk_file.name}: {e}", exc_info=True)

    if track_count == 0:
        raise RuntimeError("No valid tracks created – nothing to plot")

    logging.info(f"Created {track_count} tracks")
    logging.info("Rendering HTML")

    gv.savefig_html(out_html)
    logging.info(f"Saved: {out_html}")


if __name__ == "__main__":
    main()
