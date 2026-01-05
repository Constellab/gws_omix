#!/usr/bin/env python3

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

# Keep ONLY these essential columns in the final TSV
RUNINFO_HEADER = [
    "study_accession",
    "sample_accession",
    "experiment_accession",
    "run_accession",
    "scientific_name",
    "library_layout",
    "fastq_md5",
    "fastq_ftp",
    "base_count",
    "read_count",
]

class FastqDownloader:
    def __init__(self, accession: str, provider: str, cpus: int | str,
                 out_fastq_dir: str | Path, summary_dir: str | Path) -> None:
        self.accession = str(accession).strip()
        self.provider  = str(provider or "ena").strip().lower()
        if self.provider not in ("ena", "sra"):
            self.provider = "ena"
        self.cpus = self._to_pos_int(cpus, 1)

        self.out_fastq_dir = Path(out_fastq_dir).expanduser().resolve()
        self.summary_dir   = Path(summary_dir).expanduser().resolve()
        self.summary_tsv   = self.summary_dir / "fastq-run-info.tsv"

        self.out_fastq_dir.mkdir(parents=True, exist_ok=True)
        self.summary_dir.mkdir(parents=True, exist_ok=True)
        self._ensure_header_tsv(self.summary_tsv)

    # ---------- helpers ----------

    @staticmethod
    def _to_pos_int(value: int | str, default: int = 1) -> int:
        try:
            n = int(value)
            return n if n > 0 else default
        except Exception:
            return default

    @staticmethod
    def _ensure_header_tsv(path: Path) -> None:
        """Guarantee a header-only TSV exists."""
        if path.exists() and path.stat().st_size > 0:
            return
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\t".join(RUNINFO_HEADER) + "\n", encoding="utf-8")

    @staticmethod
    def _find_any_runinfo(out_dir: Path) -> Path | None:
        """Find any '*-run-info.tsv' produced by fastq-dl in out_dir."""
        candidates = sorted(out_dir.glob("*-run-info.tsv"))
        return candidates[0] if candidates else None

    # ---------- core ----------

    def build_cmd(self) -> list[str]:
        return [
            "fastq-dl",
            "--accession", self.accession,
            "--provider", self.provider,
            "--cpus", str(self.cpus),
            "--outdir", str(self.out_fastq_dir),
            "--prefix", "fastq",
        ]

    def run(self) -> int:
        print(f"[FASTQ-DL] accession={self.accession} provider={self.provider} cpus={self.cpus}")
        print(f"[FASTQ-DL] FASTQ out dir={self.out_fastq_dir}")
        print(f"[FASTQ-DL] Summary dir  ={self.summary_dir}")

        # 1) Check CLI presence
        if shutil.which("fastq-dl") is None:
            print("[ERROR] `fastq-dl` not found on PATH. Please install it via Bioconda.", flush=True)
            self._final_print(exit_code=127)
            return 127

        # 2) Run fastq-dl
        cmd = self.build_cmd()
        print("[FASTQ-DL] Running:", " ".join(cmd), flush=True)
        proc = subprocess.run(cmd, check=False, cwd=self.out_fastq_dir, capture_output=True, text=True)
        if proc.stdout:
            print(proc.stdout, end="")
        if proc.stderr:
            print(proc.stderr, end="")

        # 3) Move/copy produced run-info into summary dir (and remove from FASTQ dir)
        produced_tsv = self._find_any_runinfo(self.out_fastq_dir)
        if produced_tsv and produced_tsv.exists():
            try:
                if self.summary_tsv.exists():
                    self.summary_tsv.unlink()
                try:
                    produced_tsv.replace(self.summary_tsv)  # move
                except Exception:
                    self.summary_tsv.write_bytes(produced_tsv.read_bytes())  # copy
                    try:
                        produced_tsv.unlink()
                    except Exception:
                        pass
            except Exception as e:
                print(f"[WARN] Could not move/copy run-info TSV: {e}")
                # fall through; header-only already present if move/copy failed

        # 4) Ensure the summary TSV exists, then shrink it to essentials
        self._ensure_header_tsv(self.summary_tsv)
        self._shrink_runinfo_columns(self.summary_tsv)

        # 5) Final paths
        self._final_print(exit_code=proc.returncode)
        return proc.returncode

    def _shrink_runinfo_columns(self, tsv_path: Path) -> None:
        """
        Read the TSV at tsv_path and rewrite it keeping ONLY RUNINFO_HEADER columns,
        in that exact order. Missing columns are created empty. Extra columns are dropped.
        """
        try:
            text = tsv_path.read_text(encoding="utf-8")
        except Exception as e:
            print(f"[WARN] Could not read TSV for shrinking: {e}")
            self._ensure_header_tsv(tsv_path)
            return

        lines = [ln.rstrip("\n") for ln in text.splitlines()]
        if not lines:
            self._ensure_header_tsv(tsv_path)
            return

        src_header = [h.strip() for h in lines[0].split("\t")] if lines else []
        col_idx = {h: i for i, h in enumerate(src_header)}
        # Build new content
        out_lines: list[str] = []
        out_lines.append("\t".join(RUNINFO_HEADER))
        for ln in lines[1:]:
            if not ln.strip():
                continue
            parts = ln.split("\t")
            row_vals = []
            for h in RUNINFO_HEADER:
                i = col_idx.get(h)
                row_vals.append(parts[i] if i is not None and i < len(parts) else "")
            out_lines.append("\t".join(row_vals))

        try:
            tsv_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
        except Exception as e:
            print(f"[WARN] Could not rewrite TSV with essential columns: {e}")

    def _final_print(self, *, exit_code: int) -> None:
        print(f"FASTQ_DIR={self.out_fastq_dir}")
        print(f"RUN_INFO_TSV={self.summary_tsv}")
        print(f"EXIT_CODE={exit_code}")


if __name__ == "__main__":
    runner = FastqDownloader(
        accession=sys.argv[1],
        provider=sys.argv[2],
        cpus=sys.argv[3],
        out_fastq_dir=sys.argv[4],
        summary_dir=sys.argv[5],
    )
    sys.exit(runner.run())
