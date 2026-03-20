#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank, Gff
from pygenomeviz.align import MUMmer, Blast
from Bio import SeqIO


def setup_logger():
    logging.basicConfig(
        level=logging.INFO,
        format="[%(levelname)s] %(message)s",
    )


def parse_args():
    parser = argparse.ArgumentParser(description="Genome visualization using pyGenomeViz")
    parser.add_argument("--in", dest="input_path", required=True)
    parser.add_argument("--out", dest="output_dir", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--compare", choices=["mummer", "blast", "none"], default="none",
                        help="Enable genome comparison (mummer, blast, or none)")
    parser.add_argument("--min-length", type=int, default=100,
                        help="Minimum alignment length for comparison links")
    return parser.parse_args()


def collect_genome_files(input_path: Path):
    """Collect GenBank or GFF files from input path (file or folder)"""
    genbank_exts = {".gb", ".gbk", ".genbank"}
    gff_exts = {".gff", ".gff3"}

    if input_path.is_file():
        file_ext = input_path.suffix.lower()
        if file_ext not in (genbank_exts | gff_exts):
            raise RuntimeError(
                f"Unsupported file extension: {input_path.name}. "
                f"Expected GenBank ({', '.join(genbank_exts)}) or GFF ({', '.join(gff_exts)})"
            )
        return [input_path]

    if input_path.is_dir():
        genome_files = []
        # Collect GenBank files
        for pattern in ("*.gb", "*.gbk", "*.genbank"):
            genome_files.extend(input_path.glob(pattern))
        # Collect GFF files
        for pattern in ("*.gff", "*.gff3"):
            genome_files.extend(input_path.glob(pattern))

        genome_files = sorted(set(genome_files))
        if not genome_files:
            raise RuntimeError("No GenBank or GFF files found in directory")
        return genome_files

    raise RuntimeError(f"Input path does not exist: {input_path}")


def main():
    setup_logger()
    args = parse_args()

    input_path = Path(args.input_path)
    output_dir = Path(args.output_dir)
    prefix = args.prefix

    output_dir.mkdir(parents=True, exist_ok=True)
    out_html = output_dir / f"{prefix}.html"

    logging.info(f"Input path: {input_path}")
    logging.info(f"Output directory: {output_dir}")
    logging.info(f"Prefix: {prefix}")

    genome_files = collect_genome_files(input_path)

    logging.info(f"Found {len(genome_files)} genome file(s)")

    # Create GenomeViz with interactive zoom enabled
    gv = GenomeViz()
    gv.set_scale_bar(ymargin=0.5)

    track_count = 0
    genome_data = []  # Store (genome_file, track) for comparisons

    for genome_file in genome_files:
        try:
            file_ext = genome_file.suffix.lower()

            # Parse based on file extension
            if file_ext in {".gb", ".gbk", ".genbank"}:
                parser = Genbank(genome_file)
                file_type = "GenBank"
            elif file_ext in {".gff", ".gff3"}:
                parser = Gff(genome_file)
                file_type = "GFF"
            else:
                logging.warning(f"Skipping unsupported file: {genome_file.name}")
                continue

            logging.info(
                f"Loaded {file_type} file: {genome_file.name}, genome length: {parser.genome_length}"
            )

            track = gv.add_feature_track(
                parser.name or genome_file.stem,
                parser.genome_length,
            )

            track.add_sublabel()

            features = parser.extract_features(feature_type="CDS")

            if not features:
                logging.warning(f"No CDS features in {genome_file.name}")
                continue

            track.add_features(features)
            genome_data.append((genome_file, track))
            track_count += 1

        except Exception as e:
            logging.warning(f"Failed to process {genome_file.name}: {e}", exc_info=True)

    if track_count == 0:
        raise RuntimeError("No valid tracks created – nothing to plot")

    logging.info(f"Created {track_count} tracks")

    # Perform genome comparisons if requested and there are multiple genomes
    if args.compare != "none" and len(genome_data) >= 2:
        logging.info(f"Performing genome comparisons using {args.compare.upper()}")

        # Compare consecutive genome pairs
        for i in range(len(genome_data) - 1):
            genome_file1, track1 = genome_data[i]
            genome_file2, track2 = genome_data[i + 1]

            try:
                logging.info(f"Comparing {genome_file1.name} <-> {genome_file2.name}")

                # Extract sequences to temporary FASTA files
                fasta1 = output_dir / f"{genome_file1.stem}.fasta"
                fasta2 = output_dir / f"{genome_file2.stem}.fasta"

                # Parse files and extract sequences
                file_ext1 = genome_file1.suffix.lower()
                if file_ext1 in {".gb", ".gbk", ".genbank"}:
                    parser1 = Genbank(genome_file1)
                else:
                    parser1 = Gff(genome_file1)

                file_ext2 = genome_file2.suffix.lower()
                if file_ext2 in {".gb", ".gbk", ".genbank"}:
                    parser2 = Genbank(genome_file2)
                else:
                    parser2 = Gff(genome_file2)

                # Write FASTA files
                logging.info(f"Extracting sequences to FASTA...")

                # Extract sequences from parser records and write to FASTA
                with open(fasta1, 'w') as f1:
                    SeqIO.write(parser1.records, f1, "fasta")

                with open(fasta2, 'w') as f2:
                    SeqIO.write(parser2.records, f2, "fasta")

                # Verify FASTA files exist
                if not fasta1.exists() or not fasta2.exists():
                    logging.error(f"Failed to create FASTA files")
                    continue

                logging.info(f"FASTA files created: {fasta1.name} ({fasta1.stat().st_size} bytes), {fasta2.name} ({fasta2.stat().st_size} bytes)")

                # Run alignment
                logging.info(f"Running {args.compare.upper()} alignment...")
                if args.compare == "mummer":
                    # MUMmer expects a list of FASTA files
                    align = MUMmer([fasta1, fasta2], seqtype="nucleotide")
                    align.run()
                elif args.compare == "blast":
                    # BLAST expects a list of FASTA files
                    align = Blast([fasta1, fasta2], seqtype="nucleotide")
                    align.run()

                logging.info(f"Alignment completed, retrieving coordinates with min_length={args.min_length}")

                # Add alignment links to visualization
                align_coords = align.get_coords(min_length=args.min_length)

                logging.info(f"Retrieved {len(align_coords) if align_coords else 0} alignment coordinates")

                if align_coords:
                    # Debug: print first few coords
                    for idx, coord in enumerate(align_coords[:3]):
                        logging.info(f"  Alignment {idx+1}: {coord}")

                    # Add links between tracks
                    for link_data in align_coords:
                        gv.add_link(link_data[0], link_data[1], v=link_data[2], color="grey", alpha=0.5)

                    logging.info(f"✓ Added {len(align_coords)} comparison links")
                else:
                    logging.warning(f"✗ No alignments found (try lowering min_length from {args.min_length})")

            except Exception as e:
                logging.error(f"Failed to compare {genome_file1.name} with {genome_file2.name}: {e}", exc_info=True)

    elif args.compare != "none" and len(genome_data) < 2:
        logging.warning("Genome comparison requires at least 2 genomes, skipping comparison")

    logging.info("Rendering interactive HTML with zoom/pan controls")

    # Save as interactive HTML (Plotly-based with built-in zoom/pan)
    gv.savefig_html(out_html)

    logging.info(f"Saved: {out_html}")
    logging.info("💡 Use mouse wheel to zoom, click+drag to pan, double-click to reset view")


if __name__ == "__main__":
    main()