#!/usr/bin/env python3
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

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
    def __init__(
        self,
        accessions: str,
        provider: str,
        cpus: int | str,
        out_fastq_dir: str | Path,
        summary_dir: str | Path,
    ) -> None:
        raw = str(accessions).strip()
        self.accessions = [a.strip() for a in raw.split(",") if a.strip()]
        if not self.accessions:
            raise ValueError("No accession provided (expected comma-separated list or single accession).")

        self.provider = str(provider or "ena").strip().lower()
        if self.provider not in ("ena", "sra"):
            self.provider = "ena"

        self.cpus = self._to_pos_int(cpus, 1)

        self.out_fastq_dir = Path(out_fastq_dir).expanduser().resolve()
        self.summary_dir = Path(summary_dir).expanduser().resolve()
        self.summary_tsv = self.summary_dir / "fastq-run-info.tsv"

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
        if path.exists() and path.stat().st_size > 0:
            return
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\t".join(RUNINFO_HEADER) + "\n", encoding="utf-8")

    @staticmethod
    def _find_runinfo_for_prefix(out_dir: Path, prefix: str) -> Path | None:
        """Prefer the run-info for the given prefix; fallback to any '*-run-info.tsv'."""
        exact = out_dir / f"{prefix}-run-info.tsv"
        if exact.exists():
            return exact
        candidates = sorted(out_dir.glob(f"{prefix}*-run-info.tsv"))
        if candidates:
            return candidates[0]
        any_candidates = sorted(out_dir.glob("*-run-info.tsv"))
        return any_candidates[0] if any_candidates else None

    @staticmethod
    def _load_existing_run_accessions(tsv_path: Path) -> set[str]:
        if not tsv_path.exists() or tsv_path.stat().st_size == 0:
            return set()
        try:
            lines = tsv_path.read_text(encoding="utf-8").splitlines()
        except Exception:
            return set()
        if not lines:
            return set()

        header = lines[0].split("\t")
        try:
            idx = header.index("run_accession")
        except ValueError:
            return set()

        out: set[str] = set()
        for ln in lines[1:]:
            if not ln.strip():
                continue
            parts = ln.split("\t")
            if idx < len(parts) and parts[idx].strip():
                out.add(parts[idx].strip())
        return out


    def build_cmd(self, accession: str, outdir: Path, prefix: str) -> list[str]:
        return [
            "fastq-dl",
            "--accession", accession,
            "--provider", self.provider,
            "--cpus", str(self.cpus),
            "--outdir", str(outdir),
            "--prefix", prefix,
        ]

    def run(self) -> int:
        print(f"[FASTQ-DL] provider={self.provider} cpus={self.cpus}")
        print(f"[FASTQ-DL] FASTQ out dir={self.out_fastq_dir}")
        print(f"[FASTQ-DL] Summary dir  ={self.summary_dir}")
        print(f"[FASTQ-DL] Accessions   ={', '.join(self.accessions)}")

        if shutil.which("fastq-dl") is None:
            print("[ERROR] `fastq-dl` not found on PATH. Please install it via Bioconda.", flush=True)
            self._final_print(exit_code=127)
            return 127

        seen_runs = self._load_existing_run_accessions(self.summary_tsv)
        overall_rc = 0

        for acc in self.accessions:
            prefix = f"fastq_{acc}"

            print(f"\n[FASTQ-DL] ---- Processing accession={acc} ----")
            cmd = self.build_cmd(accession=acc, outdir=self.out_fastq_dir, prefix=prefix)
            print("[FASTQ-DL] Running:", " ".join(cmd), flush=True)

            proc = subprocess.run(
                cmd,
                check=False,
                cwd=self.out_fastq_dir,
                capture_output=True,
                text=True,
            )
            if proc.stdout:
                print(proc.stdout, end="")
            if proc.stderr:
                print(proc.stderr, end="")

            if proc.returncode != 0 and overall_rc == 0:
                overall_rc = proc.returncode

            produced_tsv = self._find_runinfo_for_prefix(self.out_fastq_dir, prefix=prefix)
            if produced_tsv and produced_tsv.exists():
                self._shrink_runinfo_columns(produced_tsv)
                self._append_runinfo_rows(src_tsv=produced_tsv, dest_tsv=self.summary_tsv, seen_runs=seen_runs)

                # Optionnel: supprimer le TSV local par accession pour ne garder que le global
                try:
                    produced_tsv.unlink()
                except Exception:
                    pass
            else:
                print(f"[WARN] No run-info TSV found for accession={acc} in {self.out_fastq_dir}")

        self._ensure_header_tsv(self.summary_tsv)
        self._final_print(exit_code=overall_rc)
        return overall_rc

    def _append_runinfo_rows(self, src_tsv: Path, dest_tsv: Path, seen_runs: set[str]) -> None:
        self._ensure_header_tsv(dest_tsv)

        try:
            src_lines = src_tsv.read_text(encoding="utf-8").splitlines()
        except Exception as e:
            print(f"[WARN] Could not read per-accession TSV for merge: {e}")
            return

        if len(src_lines) <= 1:
            return

        run_idx = RUNINFO_HEADER.index("run_accession")

        out_rows: list[str] = []
        for ln in src_lines[1:]:
            if not ln.strip():
                continue
            parts = ln.split("\t")
            run_acc = parts[run_idx].strip() if run_idx < len(parts) else ""
            if run_acc and run_acc in seen_runs:
                continue
            if run_acc:
                seen_runs.add(run_acc)
            out_rows.append(ln)

        if not out_rows:
            return

        try:
            with dest_tsv.open("a", encoding="utf-8") as f:
                for r in out_rows:
                    f.write(r + "\n")
        except Exception as e:
            print(f"[WARN] Could not append rows into global TSV: {e}")

    def _shrink_runinfo_columns(self, tsv_path: Path) -> None:
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

        src_header = [h.strip() for h in lines[0].split("\t")]
        col_idx = {h: i for i, h in enumerate(src_header)}

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
        print(f"\nFASTQ_DIR={self.out_fastq_dir}")
        print(f"RUN_INFO_TSV={self.summary_tsv}")
        print(f"EXIT_CODE={exit_code}")


if __name__ == "__main__":
    runner = FastqDownloader(
        accessions=sys.argv[1],
        provider=sys.argv[2],
        cpus=sys.argv[3],
        out_fastq_dir=sys.argv[4],
        summary_dir=sys.argv[5],
    )
    sys.exit(runner.run())
