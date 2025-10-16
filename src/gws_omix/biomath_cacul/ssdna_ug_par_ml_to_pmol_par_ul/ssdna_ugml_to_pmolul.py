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
    "ssDNA_ugml_to_pmolul",
    human_name="ssDNA: µg/ml to pmol/µl",
    short_description="Convert ssDNA (oligo) concentration from µg/ml to pmol/µl for a given length (bases)."
)
class SSDNAUgmlToPmolul(Task):
    """
    ssDNA: µg/ml → pmol/µl

    This module converts a single-stranded DNA (oligo) concentration from micrograms per milliliter (µg/ml)
    to picomoles per microliter (pmol/µl), given the oligo length in nucleotides (bases).

    Formula:
        pmol/µl = (µg/ml × 1000) / (N × 330)

    Where:
        - N     = oligo length (number of nucleotides, bases)
        - 330   = average molecular weight of a single nucleotide (pg/pmol) for ssDNA
        - 1000  = ml → µl conversion (1 ml = 1000 µl)

    Derivation:
        (µg/ml) × (ml / 1000 µl) × (10^6 pg / 1 µg) × (1 pmol / 330 pg) × (1 / N)
        = (µg/ml) × 1000 / (330 × N)   pmol/µl

    Example:
        For an oligo of 25 bases at 10 µg/ml:
            pmol/µl = (10 × 1000) / (25 × 330) ≈ 1.212 pmol/µl
    """
    input_specs: Final[InputSpecs] = InputSpecs({})
    output_specs: Final[OutputSpecs] = OutputSpecs({
        "json": OutputSpec(File, human_name="Result JSON")
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(default_value="ssdna_ugml_to_pmolul", short_description="Output JSON prefix"),
        "length_nt": StrParam(default_value="", short_description="Oligo Length (bases)"),
        "conc_ugml": StrParam(default_value="", short_description="Concentration (µg/ml)")
    })

    def _opt(self, name: str, val: str) -> str:
        v = (val or "").strip()
        return f' --{name} "{v}"' if v else ""

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        shell: ShellProxy = DsdnaUgToPmolShellProxyHelper.create_proxy(self.message_dispatcher)

        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "ssdna_ugml_to_pmolul"
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix = (params["prefix"] or "ssdna_ugml_to_pmolul").strip() or "ssdna_ugml_to_pmolul"
        worker = Path(__file__).with_name("_ssdna_ugml_to_pmolul.py")
        if not worker.exists():
            raise FileNotFoundError(f"Worker not found: {worker}")

        cmd = f'python3 "{worker}" --out "{out_dir}" --prefix "{prefix}"'
        cmd += self._opt("length-nt", params["length_nt"])
        cmd += self._opt("conc-ugml", params["conc_ugml"])

        print("[DEBUG] worker cmd:", cmd)
        rc = shell.run(cmd, shell_mode=True)

        json_path = out_dir / f"{prefix}.json"

        # If error, still try to return the JSON error if present
        if rc and json_path.is_file():
            return {"json": File(str(json_path))}
        if rc:
            raise RuntimeError(f"ssDNA_ugml_to_pmolul worker failed (rc={rc}).")

        if not json_path.is_file():
            raise FileNotFoundError(f"JSON not produced: {json_path}")

        return {"json": File(str(json_path))}
