#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from typing import Final

from gws_core import (
    ConfigParams, ConfigSpecs, File,
    InputSpec, InputSpecs, OutputSpec, OutputSpecs,
    ShellProxy, StrParam, Task, TaskInputs, TaskOutputs,
    task_decorator,
)

from ..dsdna_ug_to_pmol.dsdna_ug_to_pmol_env import DsdnaUgToPmolShellProxyHelper


@task_decorator(
    "dsDNA_pmol_to_ug",
    human_name="dsDNA: pmol to µg",
    short_description="Convert a dsDNA amount (pmol) to mass (µg) for a fragment of given length (bp)."
)
class DSDNAPmolToUg(Task):
    """
    dsDNA: pmol → µg

    This module converts a molecular amount of double-stranded DNA (pmol)
    into a mass (µg) for a DNA fragment of a given length in base pairs (bp).

    In other words: given N base pairs and a certain amount in picomoles,
    what is the corresponding mass in micrograms?

    Formula:
        µg = (pmol × bp × 660) / 1e6

    Where:
        - 660 = average molecular weight of one base pair of dsDNA (pg/pmol)
        - 1e6 = conversion from picograms (pg) to micrograms (µg)
        - bp  = fragment length in base pairs

    Example:
        For 1 pmol of a 3000 bp fragment:
            µg = (1 × 3000 × 660) / 1e6 = 1.98 µg

    """
    input_specs: Final[InputSpecs] = InputSpecs({})
    output_specs: Final[OutputSpecs] = OutputSpecs({
        "json": OutputSpec(File, human_name="Result JSON")
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(default_value="dsdna_pmol_to_ug", short_description="Output JSON prefix"),
        "length_bp": StrParam(default_value="", short_description="DNA Length (bp)"),
        "pmol": StrParam(default_value="", short_description="DNA Amount (pmol)")
    })

    def _opt(self, name: str, val: str) -> str:
        v = (val or "").strip()
        return f' --{name} "{v}"' if v else ""

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        shell: ShellProxy = DsdnaUgToPmolShellProxyHelper.create_proxy(self.message_dispatcher)

        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "dsdna_pmol_to_ug"
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix = (params["prefix"] or "dsdna_pmol_to_ug").strip() or "dsdna_pmol_to_ug"
        worker = Path(__file__).with_name("_dsdna_pmol_to_ug.py")
        if not worker.exists():
            raise FileNotFoundError(f"Worker not found: {worker}")

        cmd = f'python3 "{worker}" --out "{out_dir}" --prefix "{prefix}"'
        cmd += self._opt("length-bp", params["length_bp"])
        cmd += self._opt("pmol", params["pmol"])

        print("[DEBUG] worker cmd:", cmd)
        rc = shell.run(cmd, shell_mode=True)

        json_path = out_dir / f"{prefix}.json"

        if rc and json_path.is_file():
            return {"json": File(str(json_path))}
        if rc:
            raise RuntimeError(f"dsDNA_pmol_to_ug worker failed (rc={rc}).")

        if not json_path.is_file():
            raise FileNotFoundError(f"JSON not produced: {json_path}")

        return {"json": File(str(json_path))}
