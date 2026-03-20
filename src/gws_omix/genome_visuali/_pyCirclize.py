#!/usr/bin/env python3

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio.Seq import UndefinedSequenceError
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from pycirclize import Circos
from pycirclize.parser import Genbank, Gff


GENBANK_EXTENSIONS = {".gb", ".gbk", ".genbank", ".gbff"}
GFF_EXTENSIONS = {".gff", ".gff3"}
SUPPORTED_EXTENSIONS = GENBANK_EXTENSIONS | GFF_EXTENSIONS

COLORS = {
    "cds_forward": "#ff4d4d",
    "cds_reverse": "#4d5cff",
    "rrna": "#7cd67c",
    "trna": "#ff66ff",
    "misc_rna": "#d08cff",
    "gc_pos": "#666666",
    "gc_neg": "#cfcfcf",
    "skew_pos": "#9a9a00",
    "skew_neg": "#8e24aa",
    "axis_bg": "#ffffff",
    "axis_edge": "#bfbfbf",
    "tick": "#555555",
    "label": "#222222",
}

MIN_FULL_GENOME_SIZE_WARNING = 50_000


def setup_logger() -> None:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Circular genome visualization using pyCirclize"
    )
    parser.add_argument("--in", dest="input_path", required=True)
    parser.add_argument("--out", dest="output_dir", required=True)
    return parser.parse_args()


def collect_input_files(input_path: Path) -> List[Path]:
    if input_path.is_file():
        if input_path.suffix.lower() not in SUPPORTED_EXTENSIONS:
            raise RuntimeError(f"Unsupported input file extension: {input_path.name}")
        return [input_path]

    if input_path.is_dir():
        files: List[Path] = []
        for pattern in ("*.gb", "*.gbk", "*.genbank", "*.gbff", "*.gff", "*.gff3"):
            files.extend(input_path.glob(pattern))
        files = sorted(set(files))
        if not files:
            raise RuntimeError("No supported GenBank/GFF3 files found")
        return files

    raise RuntimeError(f"Input path does not exist: {input_path}")


def load_parser(file_path: Path):
    ext = file_path.suffix.lower()
    if ext in GENBANK_EXTENSIONS:
        return Genbank(file_path), "genbank"
    if ext in GFF_EXTENSIONS:
        return Gff(file_path), "gff"
    raise RuntimeError(f"Unsupported file extension: {file_path.name}")


def sanitize_text(text: str) -> str:
    return str(text).strip().replace("\n", " ").replace("\t", " ")


def short_label(text: str, max_len: int = 28) -> str:
    text = sanitize_text(text)
    return text if len(text) <= max_len else text[: max_len - 3] + "..."


def choose_tick_interval(genome_size: int) -> int:
    if genome_size >= 10_000_000:
        return 1_000_000
    if genome_size >= 4_000_000:
        return 500_000
    if genome_size >= 2_000_000:
        return 250_000
    if genome_size >= 1_000_000:
        return 100_000
    if genome_size >= 300_000:
        return 50_000
    if genome_size >= 100_000:
        return 20_000
    if genome_size >= 20_000:
        return 5_000
    return 1_000


def choose_gc_window(genome_size: int) -> int:
    if genome_size >= 8_000_000:
        return 20_000
    if genome_size >= 3_000_000:
        return 10_000
    if genome_size >= 500_000:
        return 5_000
    if genome_size >= 100_000:
        return 2_000
    return 1_000


def choose_figsize(n_sectors: int, total_size: int) -> Tuple[int, int]:
    if total_size >= 8_000_000 or n_sectors >= 8:
        return (12, 12)
    if total_size >= 2_000_000 or n_sectors >= 3:
        return (10, 10)
    return (8, 8)


def _safe_seq_len(seq) -> Optional[int]:
    if seq is None:
        return None
    try:
        return len(seq)
    except Exception:
        return None


def _safe_seq_str(seq) -> Optional[str]:
    if seq is None:
        return None

    try:
        seq_str = str(seq).upper()
    except UndefinedSequenceError:
        return None
    except Exception:
        return None

    if not seq_str:
        return None
    if set(seq_str) <= {"N"}:
        return None

    valid_bases = {b for b in seq_str if b in "ACGTN"}
    if not valid_bases:
        return None

    return seq_str


def get_seqid2size(parser, parser_type: str) -> Dict[str, int]:
    try:
        seqid2size = parser.get_seqid2size()
        if seqid2size:
            return {str(k): int(v) for k, v in seqid2size.items()}
    except Exception:
        pass

    if parser_type == "gff":
        raise RuntimeError("No sequence identifiers / sizes found in GFF")

    # Fallback: manually inspect GenBank records
    records = getattr(parser, "records", None)
    if not records:
        raise RuntimeError("No records found in GenBank parser")

    seqid2size: Dict[str, int] = {}

    for rec in records:
        seqid = str(getattr(rec, "id", None) or getattr(rec, "name", "unknown"))
        size: Optional[int] = None

        seq = getattr(rec, "seq", None)
        size = _safe_seq_len(seq)

        if size is None:
            annotations = getattr(rec, "annotations", {}) or {}
            for key in ("sequence_length", "contig_length", "size", "length"):
                value = annotations.get(key)
                if value is None:
                    continue
                try:
                    size = int(value)
                    break
                except Exception:
                    pass

        if size is None:
            features = getattr(rec, "features", []) or []
            max_end = 0
            for feat in features:
                loc = getattr(feat, "location", None)
                if loc is None:
                    continue
                try:
                    end = int(loc.end)
                    if end > max_end:
                        max_end = end
                except Exception:
                    continue
            if max_end > 0:
                size = max_end

        if size is None or size <= 0:
            raise RuntimeError(
                f"Unable to determine sequence length for GenBank record '{seqid}'"
            )

        seqid2size[seqid] = int(size)

    return seqid2size


def try_get_seq_record(parser, seqid: str):
    for attr in ("records", "_records", "seqid2record", "_seqid2record"):
        obj = getattr(parser, attr, None)
        if obj is None:
            continue

        if isinstance(obj, dict) and seqid in obj:
            return obj[seqid]

        if isinstance(obj, list):
            for rec in obj:
                rid = getattr(rec, "id", None) or getattr(rec, "name", None)
                if str(rid) == str(seqid):
                    return rec

    return None


def extract_seq_string(parser, seqid: str) -> Optional[str]:
    rec = try_get_seq_record(parser, seqid)
    if rec is None:
        return None
    seq = getattr(rec, "seq", None)
    return _safe_seq_str(seq)


def get_seqid2features_gff(parser, file_name: str) -> Dict[str, Dict[str, list]]:
    feature_types = ("CDS", "rRNA", "tRNA", "tmRNA", "ncRNA")
    out: Dict[str, Dict[str, list]] = {}

    for ft in feature_types:
        try:
            seqid2features = parser.get_seqid2features(feature_type=ft)
            out[ft] = {str(k): list(v) for k, v in seqid2features.items()}
            count = sum(len(v) for v in seqid2features.values())
            if count > 0:
                logging.info(f"{file_name}: found {count} feature(s) of type '{ft}'")
        except Exception as e:
            logging.warning(f"{file_name}: could not extract feature type '{ft}': {e}")
            out[ft] = {}

    return out


def get_seqid2features_genbank(parser, file_name: str) -> Dict[str, Dict[str, list]]:
    feature_types = ("CDS", "rRNA", "tRNA", "tmRNA", "ncRNA")
    out: Dict[str, Dict[str, list]] = {}

    for ft in feature_types:
        try:
            seqid2features = parser.get_seqid2features(feature_type=ft)
            out[ft] = {str(k): list(v) for k, v in seqid2features.items()}
            count = sum(len(v) for v in seqid2features.values())
            if count > 0:
                logging.info(f"{file_name}: found {count} feature(s) of type '{ft}'")
        except Exception as e:
            logging.warning(f"{file_name}: could not extract feature type '{ft}': {e}")
            out[ft] = {}

    return out


def gc_fraction(seq: str) -> float:
    valid = [b for b in seq if b in "ACGT"]
    if not valid:
        return 0.0
    gc = sum(1 for b in valid if b in "GC")
    return gc / len(valid)


def gc_skew(seq: str) -> float:
    g = seq.count("G")
    c = seq.count("C")
    denom = g + c
    if denom == 0:
        return 0.0
    return (g - c) / denom


def sliding_gc_metrics(
    seq: str,
    window: int,
) -> Tuple[List[int], List[float], List[float]]:
    if len(seq) < window:
        window = max(100, len(seq))

    step = max(1, window // 2)
    xs: List[int] = []
    gc_vals: List[float] = []
    skew_vals: List[float] = []

    for start in range(0, len(seq), step):
        end = min(len(seq), start + window)
        chunk = seq[start:end]
        if not chunk:
            continue

        xs.append(min(len(seq), start + len(chunk) // 2))
        gc_vals.append(gc_fraction(chunk))
        skew_vals.append(gc_skew(chunk))

    gc_avg = sum(gc_vals) / len(gc_vals) if gc_vals else 0.0
    gc_delta = [v - gc_avg for v in gc_vals]
    return xs, gc_delta, skew_vals


def feature_strand(feature) -> int:
    strand = getattr(feature.location, "strand", None)
    return -1 if strand == -1 else 1


def split_by_strand(features: list) -> Tuple[list, list]:
    plus = []
    minus = []

    for feat in features:
        if feature_strand(feat) == -1:
            minus.append(feat)
        else:
            plus.append(feat)

    return plus, minus


def safe_sector_name(
    base_name: str,
    file_stem: str,
    existing: Dict[str, dict],
    force_prefix: bool,
) -> str:
    base_name = sanitize_text(base_name)
    name = f"{file_stem}:{base_name}" if force_prefix else base_name

    if name not in existing:
        return name

    idx = 2
    while f"{name}_{idx}" in existing:
        idx += 1
    return f"{name}_{idx}"


def add_ticks(track, genome_size: int) -> None:
    interval = choose_tick_interval(genome_size)
    track.xticks_by_interval(
        interval=interval,
        label_size=7,
        outer=True,
        show_bottom_line=True,
        tick_length=1.5,
        line_kws={"ec": COLORS["tick"], "lw": 0.6},
        text_kws={"color": COLORS["label"]},
        label_orientation="vertical",
        label_formatter=lambda v: (
            f"{v / 1_000_000:.1f} Mb" if genome_size >= 1_000_000 else f"{int(v):,}"
        ),
    )


def add_feature_track(track, features: list, color: str) -> None:
    if not features:
        return

    for feature in features:
        track.genomic_features(
            feature,
            plotstyle="box",
            fc=color,
            ec="none",
            lw=0,
        )


def add_feature_track_arrow(track, features: list, color: str, alpha: float = 1.0, plotstyle: str = "arrow") -> None:
    if not features:
        return

    for feature in features:
        track.genomic_features(
            feature,
            plotstyle=plotstyle,
            fc=color,
            ec="none",
            lw=0,
            alpha=alpha,
        )


def choose_density_mode(n_cds: int) -> tuple:
    """Return (plotstyle, alpha, show_labels) based on CDS count."""
    if n_cds >= 500:
        return "box", 0.70, False
    if n_cds >= 150:
        return "arrow", 0.85, False
    return "arrow", 1.0, True


def build_dataset(files: List[Path]) -> Dict[str, dict]:
    sectors: Dict[str, dict] = {}

    for file_path in files:
        parser, parser_type = load_parser(file_path)
        seqid2size = get_seqid2size(parser, parser_type)

        if parser_type == "gff":
            seqid2features = get_seqid2features_gff(parser, file_path.name)
        else:
            seqid2features = get_seqid2features_genbank(parser, file_path.name)

        multi_seq = len(seqid2size) > 1

        for seqid, size in seqid2size.items():
            sector_name = safe_sector_name(
                base_name=seqid,
                file_stem=file_path.stem,
                existing=sectors,
                force_prefix=(multi_seq or seqid in sectors),
            )

            # Build a human-readable display label
            if parser_type == "genbank":
                rec = try_get_seq_record(parser, seqid)
                if rec is not None:
                    # rec.name = LOCUS field (e.g. "NC_002483"), rec.id = VERSION
                    locus = getattr(rec, "name", "") or getattr(rec, "id", "") or seqid
                    display_label = str(locus)
                else:
                    display_label = seqid
            else:
                display_label = seqid

            cds = list(seqid2features.get("CDS", {}).get(seqid, []))
            rrna = list(seqid2features.get("rRNA", {}).get(seqid, []))
            trna = list(seqid2features.get("tRNA", {}).get(seqid, []))
            tmrna = list(seqid2features.get("tmRNA", {}).get(seqid, []))
            ncrna = list(seqid2features.get("ncRNA", {}).get(seqid, []))

            total_annot_features = (
                len(cds) + len(rrna) + len(trna) + len(tmrna) + len(ncrna)
            )

            seq_str = extract_seq_string(parser, seqid)
            if seq_str is None and parser_type == "genbank":
                logging.info(
                    f"{file_path.name}:{seqid} sequence is undefined or unavailable; "
                    "GC content/skew tracks will be skipped"
                )

            if int(size) < MIN_FULL_GENOME_SIZE_WARNING:
                logging.warning(
                    f"{file_path.name}:{seqid} is only {int(size):,} bp. "
                    "This is a small record and is not comparable to a full bacterial chromosome plot."
                )

            logging.info(
                f"{file_path.name}:{seqid} size={int(size):,} "
                f"CDS={len(cds)} rRNA={len(rrna)} tRNA={len(trna)} "
                f"tmRNA={len(tmrna)} ncRNA={len(ncrna)}"
            )

            sectors[sector_name] = {
                "size": int(size),
                "file_stem": file_path.stem,
                "parser_type": parser_type,
                "seqid": seqid,
                "display_label": display_label,
                "sequence": seq_str,
                "cds": cds,
                "rrna": rrna,
                "trna": trna,
                "tmrna": tmrna,
                "ncrna": ncrna,
                "annot_feature_count": total_annot_features,
            }

    return sectors


def get_center_title(sectors_meta: Dict[str, dict]) -> str:
    if len(sectors_meta) == 1:
        meta = next(iter(sectors_meta.values()))
        return short_label(meta["display_label"], 48)
    return f"{len(sectors_meta)} replicons"


def plot_gc_track(sector, sequence: str, size: int) -> None:
    window = choose_gc_window(size)
    x, gc_delta, skew_vals = sliding_gc_metrics(sequence, window=window)

    gc_track = sector.add_track((60, 68))
    gc_track.axis(fc="white", ec=COLORS["axis_edge"], lw=0.5)

    gc_abs = max(max(abs(v) for v in gc_delta), 0.01)
    gc_pos = [max(v, 0.0) for v in gc_delta]
    gc_neg = [min(v, 0.0) for v in gc_delta]

    gc_track.fill_between(
        x,
        gc_pos,
        0,
        vmin=-gc_abs,
        vmax=gc_abs,
        color=COLORS["gc_pos"],
        alpha=1.0,
    )
    gc_track.fill_between(
        x,
        gc_neg,
        0,
        vmin=-gc_abs,
        vmax=gc_abs,
        color=COLORS["gc_neg"],
        alpha=1.0,
    )
    gc_track.line(
        x,
        gc_delta,
        vmin=-gc_abs,
        vmax=gc_abs,
        color="black",
        lw=0.4,
    )

    skew_track = sector.add_track((50, 58))
    skew_track.axis(fc="white", ec=COLORS["axis_edge"], lw=0.5)

    skew_abs = max(max(abs(v) for v in skew_vals), 0.01)
    skew_pos = [max(v, 0.0) for v in skew_vals]
    skew_neg = [min(v, 0.0) for v in skew_vals]

    skew_track.fill_between(
        x,
        skew_pos,
        0,
        vmin=-skew_abs,
        vmax=skew_abs,
        color=COLORS["skew_pos"],
        alpha=1.0,
    )
    skew_track.fill_between(
        x,
        skew_neg,
        0,
        vmin=-skew_abs,
        vmax=skew_abs,
        color=COLORS["skew_neg"],
        alpha=1.0,
    )
    skew_track.line(
        x,
        skew_vals,
        vmin=-skew_abs,
        vmax=skew_abs,
        color="black",
        lw=0.4,
    )


def plot_sector(sector, meta: dict) -> None:
    size = meta["size"]
    sequence = meta["sequence"]

    cds_plus, cds_minus = split_by_strand(meta["cds"])
    rrna_plus, rrna_minus = split_by_strand(meta["rrna"])
    trna_plus, trna_minus = split_by_strand(meta["trna"] + meta["tmrna"])
    ncrna_plus, ncrna_minus = split_by_strand(meta["ncrna"])

    outer_tick = sector.add_track((98, 100))
    outer_tick.axis(fc="white", ec=COLORS["axis_edge"], lw=0.8)
    add_ticks(outer_tick, size)

    plus_track = sector.add_track((90, 97))
    plus_track.axis(fc=COLORS["axis_bg"], ec=COLORS["axis_edge"], lw=0.5)

    minus_track = sector.add_track((82, 89))
    minus_track.axis(fc=COLORS["axis_bg"], ec=COLORS["axis_edge"], lw=0.5)

    rna_track = sector.add_track((72, 78))
    rna_track.axis(fc=COLORS["axis_bg"], ec=COLORS["axis_edge"], lw=0.5)

    add_feature_track(plus_track, cds_plus, COLORS["cds_forward"])
    add_feature_track(minus_track, cds_minus, COLORS["cds_reverse"])
    add_feature_track(rna_track, rrna_plus + rrna_minus, COLORS["rrna"])
    add_feature_track(rna_track, trna_plus + trna_minus, COLORS["trna"])
    add_feature_track(rna_track, ncrna_plus + ncrna_minus, COLORS["misc_rna"])

    if sequence:
        plot_gc_track(sector, sequence, size)


def plot_sector_arrow(sector, meta: dict, n_sectors: int = 1) -> None:
    size = meta["size"]
    sequence = meta["sequence"]

    cds_plus, cds_minus = split_by_strand(meta["cds"])
    rrna_plus, rrna_minus = split_by_strand(meta["rrna"])
    trna_plus, trna_minus = split_by_strand(meta["trna"] + meta["tmrna"])
    ncrna_plus, ncrna_minus = split_by_strand(meta["ncrna"])

    n_cds = len(meta["cds"])
    plotstyle, alpha, show_labels = choose_density_mode(n_cds)

    outer_tick = sector.add_track((98, 100))
    outer_tick.axis(fc="white", ec=COLORS["axis_edge"], lw=0.8)
    add_ticks(outer_tick, size)

    plus_track = sector.add_track((90, 97))
    plus_track.axis(fc=COLORS["axis_bg"], ec=COLORS["axis_edge"], lw=0.5)

    minus_track = sector.add_track((82, 89))
    minus_track.axis(fc=COLORS["axis_bg"], ec=COLORS["axis_edge"], lw=0.5)

    rna_track = sector.add_track((72, 78))
    rna_track.axis(fc=COLORS["axis_bg"], ec=COLORS["axis_edge"], lw=0.5)

    add_feature_track_arrow(plus_track, cds_plus, COLORS["cds_forward"], alpha=alpha, plotstyle=plotstyle)
    add_feature_track_arrow(minus_track, cds_minus, COLORS["cds_reverse"], alpha=alpha, plotstyle=plotstyle)
    add_feature_track_arrow(rna_track, rrna_plus + rrna_minus, COLORS["rrna"], alpha=alpha, plotstyle="box")
    add_feature_track_arrow(rna_track, trna_plus + trna_minus, COLORS["trna"], alpha=alpha, plotstyle="box")
    add_feature_track_arrow(rna_track, ncrna_plus + ncrna_minus, COLORS["misc_rna"], alpha=alpha, plotstyle="box")

    if show_labels:
        for feature in meta["cds"]:
            qualifiers = getattr(feature, "qualifiers", {}) or {}
            gene_name = qualifiers.get("gene", [None])[0]
            if gene_name:
                start = int(feature.location.start)
                end = int(feature.location.end)
                pos = (start + end) / 2
                try:
                    plus_track.annotate(pos, short_label(gene_name, 12), label_size=5)
                except Exception:
                    pass

    if sequence:
        plot_gc_track(sector, sequence, size)

    if n_sectors > 1:
        sector.text(short_label(meta["display_label"], 24), r=104, size=9, adjust_rotation=True)


def add_center_title_and_legend(fig, sectors_meta: Dict[str, dict]) -> None:
    ax = fig.axes[0]
    has_sequence = any(meta["sequence"] for meta in sectors_meta.values())

    title = get_center_title(sectors_meta)
    ax.text(
        0.5,
        0.60,
        title,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=11,
        color="black",
    )

    legend_handles = [
        Patch(facecolor=COLORS["cds_forward"], edgecolor="none", label="Forward CDS"),
        Patch(facecolor=COLORS["cds_reverse"], edgecolor="none", label="Reverse CDS"),
        Patch(facecolor=COLORS["rrna"], edgecolor="none", label="rRNA"),
        Patch(facecolor=COLORS["trna"], edgecolor="none", label="tRNA"),
    ]

    if has_sequence:
        legend_handles.extend(
            [
                Line2D(
                    [0],
                    [0],
                    color=COLORS["gc_pos"],
                    lw=4,
                    label="Positive GC Content",
                ),
                Line2D(
                    [0],
                    [0],
                    color=COLORS["gc_neg"],
                    lw=4,
                    label="Negative GC Content",
                ),
                Line2D(
                    [0],
                    [0],
                    color=COLORS["skew_pos"],
                    lw=4,
                    label="Positive GC Skew",
                ),
                Line2D(
                    [0],
                    [0],
                    color=COLORS["skew_neg"],
                    lw=4,
                    label="Negative GC Skew",
                ),
            ]
        )

    ax.legend(
        handles=legend_handles,
        loc="center",
        bbox_to_anchor=(0.5, 0.43),
        frameon=False,
        fontsize=9,
        ncol=1,
        handlelength=1.2,
        handletextpad=0.6,
        borderpad=0.2,
        labelspacing=0.5,
    )


def render_arrow_plot(sectors_meta: Dict[str, dict], out_png: Path) -> None:
    sectors = {name: meta["size"] for name, meta in sectors_meta.items()}
    total_size = sum(sectors.values())
    n_sectors = len(sectors)

    circos = Circos(sectors=sectors, space=0 if n_sectors == 1 else 4)

    for sector in circos.sectors:
        meta = sectors_meta[sector.name]
        plot_sector_arrow(sector, meta, n_sectors=n_sectors)

    figsize = choose_figsize(n_sectors, total_size)
    fig = circos.plotfig(dpi=300, figsize=figsize)
    add_center_title_and_legend(fig, sectors_meta)
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    logging.info(f"Saved: {out_png}")


def main() -> None:
    setup_logger()
    args = parse_args()

    input_path = Path(args.input_path)
    output_dir = Path(args.output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    files = collect_input_files(input_path)
    logging.info(f"Found {len(files)} supported input file(s)")

    skipped = 0
    for file_path in files:
        logging.info(f"Processing: {file_path.name}")
        try:
            sectors_meta = build_dataset([file_path])
        except Exception as e:
            logging.warning(f"Skipping {file_path.name}: {e}")
            skipped += 1
            continue

        if not sectors_meta:
            logging.warning(f"Skipping {file_path.name}: no valid sectors")
            skipped += 1
            continue

        total_annot_features = sum(
            meta["annot_feature_count"] for meta in sectors_meta.values()
        )
        if total_annot_features == 0:
            logging.warning(
                f"Skipping {file_path.name}: no CDS/rRNA/tRNA/tmRNA/ncRNA annotations found. "
                "A GenBank containing only 'source' and 'CONTIG' features is insufficient."
            )
            skipped += 1
            continue

        out_png = output_dir / f"{file_path.stem}.png"
        render_arrow_plot(sectors_meta, out_png)

    if skipped == len(files):
        raise RuntimeError(
            f"All {len(files)} input file(s) were skipped. "
            "No supported annotation features (CDS/rRNA/tRNA) were found in any file."
        )


if __name__ == "__main__":
    try:
        main()
    except RuntimeError as e:
        logging.shutdown()
        print(f"\nERROR: {e}", file=sys.stderr, flush=True)
        sys.exit(1)