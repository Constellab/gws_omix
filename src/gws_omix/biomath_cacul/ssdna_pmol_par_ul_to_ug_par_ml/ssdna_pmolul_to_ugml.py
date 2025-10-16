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
    "ssDNA_pmolul_to_ugml",
    human_name="ssDNA: pmol/µl to µg/ml",
    short_description="Convert ssDNA (oligo) concentration from pmol/µl to µg/ml for a given length (bases)."
)
class SSDNAPmolulToUgml(Task):

    """
    ssDNA: pmol/µl → µg/ml

    This module converts a single-stranded DNA (oligo) concentration from picomoles per microliter (pmol/µl)
    to micrograms per milliliter (µg/ml), given the oligo length in nucleotides (bases).

    Formula:
        µg/ml = (pmol/µl × 1000 µl/ml × 330 pg/pmol × N) × (1 µg / 1e6 pg)
            = (pmol/µl × N × 330) / 1000

    Where:
        - N     = oligo length (number of nucleotides, bases)
        - 330   = average molecular weight of a single nucleotide (pg/pmol) for ssDNA
        - 1000  = µl → ml conversion (1000 µl in 1 ml)
        - 1e6   = pg → µg conversion (1 µg = 10^6 pg)

    Example:
        For a 25-base oligo at 1.212 pmol/µl:
            µg/ml = (1.212 × 25 × 330) / 1000 ≈ 10.0 µg/ml
    """
    input_specs: Final[InputSpecs] = InputSpecs({})
    output_specs: Final[OutputSpecs] = OutputSpecs({
        "json": OutputSpec(File, human_name="Result JSON")
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(default_value="ssdna_pmolul_to_ugml", short_description="Output JSON prefix"),
        "length_nt": StrParam(default_value="", short_description="Oligo Length (bases)"),
        "conc_pmolul": StrParam(default_value="", short_description="Concentration (pmol/µl)")
    })

    def _opt(self, name: str, val: str) -> str:
        v = (val or "").strip()
        return f' --{name} "{v}"' if v else ""

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        shell: ShellProxy = DsdnaUgToPmolShellProxyHelper.create_proxy(self.message_dispatcher)

        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "ssdna_pmolul_to_ugml"
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix = (params["prefix"] or "ssdna_pmolul_to_ugml").strip() or "ssdna_pmolul_to_ugml"
        worker = Path(__file__).with_name("_ssdna_pmolul_to_ugml.py")
        if not worker.exists():
            raise FileNotFoundError(f"Worker not found: {worker}")

        cmd = f'python3 "{worker}" --out "{out_dir}" --prefix "{prefix}"'
        cmd += self._opt("length-nt", params["length_nt"])
        cmd += self._opt("conc-pmolul", params["conc_pmolul"])

        print("[DEBUG] worker cmd:", cmd)
        rc = shell.run(cmd, shell_mode=True)

        json_path = out_dir / f"{prefix}.json"

        # If error, still try to return the JSON error if present
        if rc and json_path.is_file():
            return {"json": File(str(json_path))}
        if rc:
            raise RuntimeError(f"ssDNA_pmolul_to_ugml worker failed (rc={rc}).")

        if not json_path.is_file():
            raise FileNotFoundError(f"JSON not produced: {json_path}")

        return {"json": File(str(json_path))}
