#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from typing import Final

from gws_core import (
    ConfigParams, ConfigSpecs, File,
    InputSpec, InputSpecs, OutputSpec, OutputSpecs,
    ShellProxy, StrParam, Task, TaskInputs, TaskOutputs,
    task_decorator
)

from ..dsdna_ug_to_pmol.dsdna_ug_to_pmol_env import DsdnaUgToPmolShellProxyHelper

@task_decorator(
    "LinearDNA_ug_to_pmol_ends",
    human_name="Linear DNA: µg to pmol of Ends",
    short_description="Convert mass (µg) of linear dsDNA to pmol of DNA ends, given length (kb)."
)
class LinearDNAUgToPmolEnds(Task):

    """
    Linear DNA: µg → pmol of Ends

    This tool converts the mass of a **linear double-stranded DNA** fragment (in micrograms, µg)
    to the **picomoles of DNA ends** available, given the fragment length in kilobases (kb).

    Rationale:
    - Molecular weight per base pair (bp) of dsDNA ≈ 660 pg/pmol
    - Length in kb → bp via: bp = kb × 1000
    - Picomoles of DNA molecules:
        pmol_molecules = (µg × 1e6 pg/µg) / (bp × 660 pg/pmol)
    - A linear DNA molecule has **2 ends**, so:
        pmol_ends = pmol_molecules × 2

    Compact formula:
        pmol_ends = (µg × 1e6) / (kb × 1000 × 660) × 2
                = (µg × 2 × 10^6) / (kb × 660 × 1000)
    """

    input_specs: Final[InputSpecs] = InputSpecs({})
    output_specs: Final[OutputSpecs] = OutputSpecs({
        "json": OutputSpec(File, human_name="Result JSON")
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(default_value="linear_dna_ug_to_pmol_ends", short_description="Output JSON prefix"),
        "mass_ug": StrParam(default_value="", short_description="µg of DNA"),
        "length_kb": StrParam(default_value="", short_description="DNA Length (kb)")
    })

    def _opt(self, name: str, val: str) -> str:
        v = (val or "").strip()
        return f' --{name} "{v}"' if v else ""

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        shell: ShellProxy = DsdnaUgToPmolShellProxyHelper.create_proxy(self.message_dispatcher)

        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "linear_dna_ug_to_pmol_ends"
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix = (params["prefix"] or "linear_dna_ug_to_pmol_ends").strip() or "linear_dna_ug_to_pmol_ends"
        worker = Path(__file__).with_name("_linear_dna_ug_to_pmol_ends.py")
        if not worker.exists():
            raise FileNotFoundError(f"Worker not found: {worker}")

        cmd = f'python3 "{worker}" --out "{out_dir}" --prefix "{prefix}"'
        cmd += self._opt("mass-ug", params["mass_ug"])
        cmd += self._opt("length-kb", params["length_kb"])

        print("[DEBUG] worker cmd:", cmd)
        rc = shell.run(cmd, shell_mode=True)

        json_path = out_dir / f"{prefix}.json"

        # If error, still try to return the JSON error if present
        if rc and json_path.is_file():
            return {"json": File(str(json_path))}
        if rc:
            raise RuntimeError(f"linear_dna_ug_to_pmol_ends worker failed (rc={rc}).")

        if not json_path.is_file():
            raise FileNotFoundError(f"JSON not produced: {json_path}")

        return {"json": File(str(json_path))}
