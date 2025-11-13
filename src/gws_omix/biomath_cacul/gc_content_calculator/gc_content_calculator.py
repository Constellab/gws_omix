#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from typing import Final
import shlex

from gws_core import (
    ConfigParams, ConfigSpecs, File,
    InputSpec, InputSpecs, OutputSpec, OutputSpecs,
    ShellProxy, StrParam, Task, TaskInputs, TaskOutputs,
    task_decorator
)

from ..dsdna_ug_to_pmol.dsdna_ug_to_pmol_env import DsdnaUgToPmolShellProxyHelper


@task_decorator(
    "gc_content_calculator",
    human_name="GC Content Calculator",
    short_description="Compute GC% (100-nt window), CpG/GpC per 100 nt, predict CpG islands, and render an interactive HTML report with sequence highlighting."
)
class GCContentCalculator(Task):
    """
    This task reads a DNA/RNA sequence (FASTA or plain text), computes per-window
    GC% (window=100 nt, step=1), CpG and GpC counts per 100 nt, predicts CpG
    islands (len≥200, GC%≥50, Obs/Exp CpG≥0.6), and produces a self-contained
    HTML report with two Plotly subplots and a formatted sequence panel.
    Selecting a region on the GC% plot highlights the corresponding bases in the
    sequence view.

    Key features:
    - Input: FASTA or plain text; non-IUPAC characters removed; case-insensitive.
    - Type autodetect: RNA if 'U' present and no 'T'; otherwise DNA.
    - Output: <out_dir>/<prefix>.gc_content.html containing:
        1) GC% (100-nt) colored bars (top subplot)
        2) CpG/GpC per 100 nt with predicted CpG islands (bottom subplot)
        3) Full sequence panel (monospace, grouped by 10/100, 1-based coordinates)
        4) Interactive link: selection on GC% → highlight in sequence
    """
    input_specs: Final[InputSpecs] = InputSpecs({
        "input_path": InputSpec(
            File,
            human_name="Sequence File",
            short_description="FASTA or plain DNA/RNA sequence (first record used)."
        )
    })

    output_specs: Final[OutputSpecs] = OutputSpecs({
        "report_html": OutputSpec(File, human_name="GC Content HTML Report")
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(default_value="gc_content", short_description="Output file prefix (no extension)")
    })

    def _q(self, x) -> str:
        return shlex.quote(str(x))

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        shell: ShellProxy = DsdnaUgToPmolShellProxyHelper.create_proxy(self.message_dispatcher)

        in_file: File = inputs["input_path"]
        in_path = Path(in_file.path)
        if not in_path.is_file():
            raise FileNotFoundError(f"Input file not found: {in_path}")

        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "gc_content_calculator"
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix = (params["prefix"] or "gc_content").strip() or "gc_content"
        worker = Path(__file__).with_name("_gc_content_calculator.py")
        if not worker.exists():
            raise FileNotFoundError(f"Worker not found: {worker}")

        cmd = (
            f"python3 {self._q(worker)}"
            f" --in {self._q(in_path)}"
            f" --out {self._q(out_dir)}"
            f" --prefix {self._q(prefix)}"
        )

        print("[DEBUG] worker cmd:", cmd)
        rc = shell.run(cmd, shell_mode=True)

        html_path = out_dir / f"{prefix}.gc_content.html"
        if html_path.is_file():
            return {"report_html": File(str(html_path))}
        if rc:
            raise RuntimeError(f"GC content worker failed (rc={rc}).")
        raise FileNotFoundError(f"Report not produced: {html_path}")
