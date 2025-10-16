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
    "coding_capacity",
    human_name="Coding Capacity of DNA",
    short_description="Enter exactly one of: DNA Length (bp), Protein Length (aa), or Protein Size (kDa); the other two are computed."
)
class CodingCapacity(Task):
    """
    Coding Capacity of DNA

    Enter exactly ONE parameter:
      • DNA Length (bp), or
      • Protein Length (amino acids), or
      • Protein Size (kDa)

    The tool computes the other two using:
      • amino acids = DNA (bp) / 3
      • DNA (bp)    = 3 × amino acids
      • protein size (kDa) = amino acids × 0.11

    Notes:
      • 0.11 kDa is the average molecular weight per amino acid.
      • If you provide multiple inputs, the worker returns a JSON error message.
    """

    input_specs: Final[InputSpecs] = InputSpecs({})
    output_specs: Final[OutputSpecs] = OutputSpecs({
        "json": OutputSpec(File, human_name="Result JSON")
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(default_value="coding_capacity", short_description="Output JSON prefix"),
        "dna_bp": StrParam(default_value="", short_description="DNA Length (bp) — provide only ONE of dna_bp/aa/protein_kda"),
        "aa": StrParam(default_value="", short_description="Protein Length (amino acids) — provide only ONE of dna_bp/aa/protein_kda"),
        "protein_kda": StrParam(default_value="", short_description="Protein Size (kDa) — provide only ONE of dna_bp/aa/protein_kda"),
    })

    def _q(self, x) -> str:
        return shlex.quote(str(x))

    def _opt(self, name: str, val: str) -> str:
        v = (val or "").strip()
        return f" --{name} {self._q(v)}" if v else ""

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        shell: ShellProxy = DsdnaUgToPmolShellProxyHelper.create_proxy(self.message_dispatcher)
        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "coding_capacity"
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix = (params["prefix"] or "coding_capacity").strip() or "coding_capacity"
        worker = Path(__file__).with_name("_coding_capacity.py")
        if not worker.exists():
            raise FileNotFoundError(f"Worker not found: {worker}")

        cmd = f"python3 {self._q(worker)} --out {self._q(out_dir)} --prefix {self._q(prefix)}"
        cmd += self._opt("dna-bp", params["dna_bp"])
        cmd += self._opt("aa", params["aa"])
        cmd += self._opt("protein-kda", params["protein_kda"])

        print("[DEBUG] worker cmd:", cmd)
        rc = shell.run(cmd, shell_mode=True)

        json_path = out_dir / f"{prefix}.json"
        # Surface JSON even if rc != 0; but worker exits 0 anyway with inline error text.
        if rc and json_path.is_file():
            return {"json": File(str(json_path))}
        if rc:
            raise RuntimeError(f"Coding Capacity worker failed (rc={rc}).")
        if not json_path.is_file():
            raise FileNotFoundError(f"JSON not produced: {json_path}")

        return {"json": File(str(json_path))}
