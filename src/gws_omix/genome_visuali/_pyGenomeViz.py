#!/usr/bin/env python3

import argparse
import logging
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Tuple

from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank, Gff

GENBANK_EXTENSIONS = {".gb", ".gbk", ".genbank", ".gbff"}
GFF_EXTENSIONS = {".gff", ".gff3"}
SUPPORTED_EXTENSIONS = GENBANK_EXTENSIONS | GFF_EXTENSIONS

COMPARISON_MODES = {
    "visualization_only",
    "genome_comparison",
    "cds_protein_homology",
}

FEATURE_STYLE_MAP: Dict[str, Dict[str, Any]] = {
    "CDS": {"fc": "#ef5350", "plotstyle": "bigarrow", "lw": 0.35},
    "gene": {"fc": "#42a5f5", "plotstyle": "arrow", "lw": 0.25},
    "mRNA": {"fc": "#66bb6a", "plotstyle": "arrow", "lw": 0.25},
    "transcript": {"fc": "#66bb6a", "plotstyle": "arrow", "lw": 0.25},
    "exon": {"fc": "#26a69a", "plotstyle": "box", "lw": 0.25},
    "rRNA": {"fc": "#ffa726", "plotstyle": "bigbox", "lw": 0.25},
    "tRNA": {"fc": "#ab47bc", "plotstyle": "bigbox", "lw": 0.25},
    "tmRNA": {"fc": "#8d6e63", "plotstyle": "box", "lw": 0.25},
    "ncRNA": {"fc": "#d4e157", "plotstyle": "box", "lw": 0.25},
}

FEATURE_TYPES = ["CDS", "gene", "mRNA", "exon", "rRNA", "tRNA", "tmRNA", "ncRNA"]


def setup_logger() -> None:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Genome visualization / comparison using pyGenomeViz")
    parser.add_argument("--in", dest="input_path", required=True)
    parser.add_argument("--out", dest="output_dir", required=True)
    parser.add_argument("--prefix", default="pygenomeviz")
    parser.add_argument("--comparison-mode", required=True, choices=sorted(COMPARISON_MODES))
    return parser.parse_args()


def collect_input_files(input_path: Path) -> List[Path]:
    if not input_path.exists():
        raise RuntimeError(f"Input path does not exist: {input_path}")

    if input_path.is_file():
        ext = input_path.suffix.lower()
        if ext not in SUPPORTED_EXTENSIONS:
            raise RuntimeError(
                f"Unsupported input file extension: {input_path.name}. "
                f"Expected one of: {', '.join(sorted(SUPPORTED_EXTENSIONS))}"
            )
        return [input_path]

    files: List[Path] = []
    for pattern in ("*.gb", "*.gbk", "*.genbank", "*.gbff", "*.gff", "*.gff3"):
        files.extend(input_path.glob(pattern))
    files = sorted(set(files))

    if not files:
        raise RuntimeError("No supported GenBank/GFF files found in the input folder.")
    return files


def split_files_by_type(files: List[Path]) -> Tuple[List[Path], List[Path]]:
    genbank_files = [f for f in files if f.suffix.lower() in GENBANK_EXTENSIONS]
    gff_files = [f for f in files if f.suffix.lower() in GFF_EXTENSIONS]
    return genbank_files, gff_files


def ensure_executable_exists(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(
            f"Required executable not found in PATH: {name}. "
            f"Please install it in the environment first."
        )


def run_command(cmd: List[str], cwd: Path) -> None:
    logging.info("Running: %s", " ".join(cmd))
    proc = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)

    if proc.stdout:
        logging.info(proc.stdout.strip())
    if proc.stderr:
        logging.warning(proc.stderr.strip())

    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {proc.returncode}:\n"
            f"{' '.join(cmd)}\n"
            f"STDERR:\n{proc.stderr}"
        )


def load_parser(file_path: Path):
    ext = file_path.suffix.lower()
    if ext in GENBANK_EXTENSIONS:
        return Genbank(file_path), "genbank"
    if ext in GFF_EXTENSIONS:
        return Gff(file_path), "gff"
    raise RuntimeError(f"Unsupported input file: {file_path.name}")


def parser_name(parser, fallback: str) -> str:
    value = getattr(parser, "name", None)
    return str(value) if value else fallback


def parser_seqid2size(parser) -> Dict[str, int]:
    if hasattr(parser, "get_seqid2size"):
        value = parser.get_seqid2size()
        if isinstance(value, dict) and value:
            return {str(k): int(v) for k, v in value.items()}
    if hasattr(parser, "genome_length"):
        return {parser_name(parser, "genome"): int(getattr(parser, "genome_length"))}
    raise RuntimeError("Unable to determine sequence sizes from parser")


def feature_sort_key(feature) -> Tuple[int, int]:
    loc = getattr(feature, "location", None)
    if loc is None:
        return (0, 0)
    return (int(loc.start), int(loc.end))


def extract_features_safe(parser, feature_type: str, seq_size: int):
    try:
        return list(parser.extract_features(feature_type=feature_type, target_range=(0, seq_size)))
    except TypeError:
        try:
            return list(parser.extract_features(feature_type, target_range=(0, seq_size)))
        except Exception:
            return []
    except Exception:
        return []


def add_feature_batch_compat(target, features, plotstyle: str, fc: str, lw: float) -> bool:
    if not features:
        return True

    for kwargs in (
        {"plotstyle": plotstyle, "fc": fc, "lw": lw},
        {"plotstyle": plotstyle, "fc": fc},
        {"plotstyle": plotstyle},
        {},
    ):
        try:
            target.add_features(features, **kwargs)
            return True
        except Exception:
            continue
    return False


def add_feature_single_compat(target, start: int, end: int, strand: int, plotstyle: str, fc: str, lw: float) -> None:
    for kwargs in (
        {"plotstyle": plotstyle, "fc": fc, "lw": lw},
        {"plotstyle": plotstyle, "fc": fc},
        {"plotstyle": plotstyle},
        {},
    ):
        try:
            target.add_feature(start, end, strand, **kwargs)
            return
        except Exception:
            continue
    target.add_feature(start, end, strand)


def add_exon_feature_compat(target, locs, strand: int, fc: str, lw: float) -> None:
    try:
        target.add_exon_feature(locs, strand=strand, plotstyle="box", fc=fc, lw=lw)
        return
    except Exception:
        pass
    try:
        target.add_exon_feature(locs, strand=strand, fc=fc, lw=lw)
        return
    except Exception:
        pass
    target.add_exon_feature(locs, strand=strand)


def build_visualization_only_html(files: List[Path], output_html: Path) -> None:
    gv = GenomeViz(
        track_align_type="left",
        fig_width=50,
        fig_track_height=max(1.0, min(2.5, 6.0 / max(1, len(files))))
    )

    for file_path in files:
        parser, _file_format = load_parser(file_path)
        seqid2size = parser_seqid2size(parser)
        track_name = file_path.stem

        multi_segment = len(seqid2size) > 1
        segments = seqid2size if multi_segment else next(iter(seqid2size.values()))
        track = gv.add_feature_track(track_name, segments, align_label=False)

        try:
            track.add_sublabel()
        except Exception:
            pass

        for seqid, seq_size in seqid2size.items():
            target = track
            if multi_segment and hasattr(track, "get_segment"):
                try:
                    target = track.get_segment(seqid)
                except Exception:
                    target = track

            for feature_type in FEATURE_TYPES:
                if feature_type == "exon":
                    continue

                features = sorted(
                    extract_features_safe(parser, feature_type, seq_size),
                    key=feature_sort_key,
                )
                if not features:
                    continue

                style = dict(FEATURE_STYLE_MAP.get(feature_type, {"fc": "#90caf9", "plotstyle": "arrow", "lw": 0.2}))
                fc = style.get("fc", "#90caf9")
                plotstyle = style.get("plotstyle", "arrow")
                lw = style.get("lw", 0.3)

                ok = add_feature_batch_compat(
                    target=target,
                    features=features,
                    plotstyle=plotstyle,
                    fc=fc,
                    lw=lw,
                )
                if not ok:
                    for feature in features:
                        loc = getattr(feature, "location", None)
                        if loc is None:
                            continue
                        start, end = int(loc.start), int(loc.end)
                        strand = -1 if getattr(loc, "strand", None) == -1 else 1
                        add_feature_single_compat(
                            target=target,
                            start=start,
                            end=end,
                            strand=strand,
                            plotstyle=plotstyle,
                            fc=fc,
                            lw=lw,
                        )

            exon_features = sorted(
                extract_features_safe(parser, "exon", seq_size),
                key=feature_sort_key,
            )
            grouped = defaultdict(list)
            for exon in exon_features:
                grouped["default"].append(exon)

            for _, exons in grouped.items():
                if not exons:
                    continue
                locs = []
                strand = 1
                for exon in exons:
                    loc = getattr(exon, "location", None)
                    if loc is None:
                        continue
                    strand = -1 if getattr(loc, "strand", None) == -1 else 1
                    parts = getattr(loc, "parts", None)
                    if parts:
                        for part in parts:
                            locs.append((int(part.start), int(part.end)))
                    else:
                        locs.append((int(loc.start), int(loc.end)))

                if locs:
                    try:
                        add_exon_feature_compat(
                            target=target,
                            locs=locs,
                            strand=strand,
                            fc=FEATURE_STYLE_MAP["exon"]["fc"],
                            lw=0.25,
                        )
                    except Exception:
                        pass

    gv.savefig_html(output_html)
    postprocess_html_width(output_html)
    hide_pseudo_column(output_html)


def build_comparison_html_with_mummer(files: List[Path], output_dir: Path, final_html: Path) -> None:
    ensure_executable_exists("pgv-mummer")
    ensure_executable_exists("nucmer")

    run_dir = output_dir / "mummer_run"
    run_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "pgv-mummer",
        *[str(f) for f in files],
        "-o", str(run_dir),
        "--formats", "html",
        "--seqtype", "nucleotide",
        "--track_align_type", "left",
        "--show_scale_xticks",
        "--curve",
        "--feature_type2color", "CDS:orange", "rRNA:lime", "tRNA:magenta",
        "--feature_linewidth", "0.3",
    ]
    run_command(cmd, cwd=output_dir)

    generated_html = find_generated_html(run_dir)
    shutil.copy2(generated_html, final_html)
    postprocess_html_width(final_html)
    hide_pseudo_column(final_html)


def build_comparison_html_with_mmseqs(files: List[Path], output_dir: Path, final_html: Path) -> None:
    ensure_executable_exists("pgv-mmseqs")
    ensure_executable_exists("mmseqs")

    run_dir = output_dir / "mmseqs_run"
    run_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "pgv-mmseqs",
        *[str(f) for f in files],
        "-o", str(run_dir),
        "--formats", "html",
        "--track_align_type", "left",
        "--show_scale_xticks",
        "--curve",
        "--feature_type2color", "CDS:skyblue",
        "--feature_linewidth", "0.3",
        "--normal_link_color", "chocolate",
        "--inverted_link_color", "limegreen",
    ]
    run_command(cmd, cwd=output_dir)

    generated_html = find_generated_html(run_dir)
    shutil.copy2(generated_html, final_html)
    postprocess_html_width(final_html)
    hide_pseudo_column(final_html)


def find_generated_html(run_dir: Path) -> Path:
    html_files = sorted(run_dir.rglob("*.html"))
    if not html_files:
        raise RuntimeError(f"No HTML file was generated in {run_dir}")
    return html_files[0]


def postprocess_html_width(html_path: Path) -> None:
    content = html_path.read_text(encoding="utf-8")
    inject = """
<style>
html, body { margin: 0; padding: 0; width: 100%; }
.bk-root { width: 100% !important; }
.bk-root > div { width: 100% !important; }
.bk { width: 100% !important; }
</style>
<script>
(function () {
  function stretchBokeh() {
    var root = document.querySelector('.bk-root');
    if (!root) return;
    root.style.width = '100%';
    window.dispatchEvent(new Event('resize'));
  }
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', function () { setTimeout(stretchBokeh, 200); });
  } else {
    setTimeout(stretchBokeh, 200);
  }
  window.addEventListener('load', function () { setTimeout(stretchBokeh, 400); });
}());
</script>
"""
    content = content.replace("</head>", inject + "</head>", 1)
    html_path.write_text(content, encoding="utf-8")


def hide_pseudo_column(html_path: Path) -> None:
    content = html_path.read_text(encoding="utf-8")
    inject = """
<script>
(function () {
  function hidePseudoColumn() {
    const tables = Array.from(document.querySelectorAll("table"));
    tables.forEach((table) => {
      const headerCells = Array.from(table.querySelectorAll("thead th, tr th"));
      if (!headerCells.length) return;

      let pseudoIndex = -1;
      headerCells.forEach((th, idx) => {
        const txt = (th.textContent || "").trim().toLowerCase();
        if (txt === "pseudo") pseudoIndex = idx;
      });

      if (pseudoIndex < 0) return;

      headerCells.forEach((th, idx) => {
        if (idx === pseudoIndex) th.style.display = "none";
      });

      const rows = table.querySelectorAll("tr");
      rows.forEach((tr) => {
        const cells = Array.from(tr.children);
        if (cells[pseudoIndex]) {
          cells[pseudoIndex].style.display = "none";
        }
      });
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", function () {
      setTimeout(hidePseudoColumn, 300);
      setTimeout(hidePseudoColumn, 1000);
    });
  } else {
    setTimeout(hidePseudoColumn, 300);
    setTimeout(hidePseudoColumn, 1000);
  }

  window.addEventListener("load", function () {
    setTimeout(hidePseudoColumn, 1200);
  });
}());
</script>
"""
    content = content.replace("</body>", inject + "\n</body>", 1)
    html_path.write_text(content, encoding="utf-8")


def main() -> None:
    setup_logger()
    args = parse_args()

    input_path = Path(args.input_path)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    prefix = (args.prefix or "pygenomeviz").strip() or "pygenomeviz"
    requested_mode = args.comparison_mode
    final_html = output_dir / f"{prefix}.html"

    files = collect_input_files(input_path)
    genbank_files, gff_files = split_files_by_type(files)

    effective_mode = requested_mode
    if gff_files:
        logging.info(
            "Detected GFF/GFF3 input. Forcing comparison mode to visualization_only."
        )
        effective_mode = "visualization_only"

    logging.info("Collected %d file(s)", len(files))
    logging.info("GenBank files: %d", len(genbank_files))
    logging.info("GFF files: %d", len(gff_files))
    logging.info("Requested mode: %s", requested_mode)
    logging.info("Effective mode: %s", effective_mode)

    if effective_mode == "visualization_only":
        build_visualization_only_html(files, final_html)
        logging.info("Visualization HTML saved to: %s", final_html)
        return

    if effective_mode == "genome_comparison":
        if len(genbank_files) < 2:
            raise RuntimeError(
                "genome_comparison requires at least 2 GenBank files."
            )
        build_comparison_html_with_mummer(genbank_files, output_dir, final_html)
        logging.info("Genome comparison HTML saved to: %s", final_html)
        return

    if effective_mode == "cds_protein_homology":
        if len(genbank_files) < 2:
            raise RuntimeError(
                "cds_protein_homology requires at least 2 GenBank files."
            )
        build_comparison_html_with_mmseqs(genbank_files, output_dir, final_html)
        logging.info("CDS / protein homology HTML saved to: %s", final_html)
        return

    raise RuntimeError(f"Unsupported comparison mode: {effective_mode}")


if __name__ == "__main__":
    try:
        main()
    except RuntimeError as e:
        logging.shutdown()
        print(f"\\nERROR: {e}", file=sys.stderr, flush=True)
        sys.exit(1)