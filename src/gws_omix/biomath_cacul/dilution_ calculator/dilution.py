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
    "dilution_calculator",
    human_name="Dilution Calculator",
    short_description="This calculator will determine the volume of stock solution required to make a working solution of a specific concentration and volume."
)
class DilutionCalculator(Task):
    """

   This calculator will determine the volume of stock solution required to make a working solution of a specific concentration and volume.

   Enter the concentration of your stock solution, and the concentration and volume of the desired working solution.


    C1 · V1 = C2 · V2   ⇒   V1 = (C2 · V2) / C1

    """

    input_specs: Final[InputSpecs] = InputSpecs({})
    output_specs: Final[OutputSpecs] = OutputSpecs({
        "json": OutputSpec(File, human_name="Result JSON")
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(default_value="dilution", short_description="Output JSON prefix"),

        # Stock
        "stock_conc": StrParam(default_value="", short_description="Stock concentration (numeric)"),
        "stock_unit": StrParam(
            default_value="mM",
            allowed_values=["M", "mM", "µM", "nM", "pM"],
            short_description="Stock concentration unit"
        ),

        # Final target
        "final_conc": StrParam(default_value="", short_description="Final concentration (numeric)"),
        "final_unit": StrParam(
            default_value="mM",
            allowed_values=["M", "mM", "µM", "nM", "pM"],
            short_description="Final concentration unit"
        ),
        "final_vol": StrParam(default_value="", short_description="Final volume (numeric)"),
        "final_vol_unit": StrParam(
            default_value="mL",
            allowed_values=["L", "mL", "µL"],
            short_description="Final volume unit"
        ),

        # NEW: output display unit for V1 (optional). If left blank, worker uses final_vol_unit.
        "stock_vol_unit_out": StrParam(
            default_value="mL",
            allowed_values=["L", "mL", "µL"],
            short_description="(Optional) Output unit for V1. Leave empty to use the same unit as final volume."
        ),
    })

    def _q(self, x) -> str:
        return shlex.quote(str(x))

    def _opt(self, name: str, val: str) -> str:
        v = (val or "").strip()
        return f" --{name} {self._q(v)}" if v else ""

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        shell: ShellProxy = DsdnaUgToPmolShellProxyHelper.create_proxy(self.message_dispatcher)

        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "dilution"
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix = (params["prefix"] or "dilution").strip() or "dilution"
        worker = Path(__file__).with_name("_dilution.py")
        if not worker.exists():
            raise FileNotFoundError(f"Worker not found: {worker}")

        cmd = f"python3 {self._q(worker)} --out {self._q(out_dir)} --prefix {self._q(prefix)}"
        cmd += self._opt("stock-conc", params["stock_conc"])
        cmd += self._opt("stock-unit", params["stock_unit"])
        cmd += self._opt("final-conc", params["final_conc"])
        cmd += self._opt("final-unit", params["final_unit"])
        cmd += self._opt("final-vol", params["final_vol"])
        cmd += self._opt("final-vol-unit", params["final_vol_unit"])
        # NEW: forward optional display unit for V1
        cmd += self._opt("stock-vol-unit-out", params["stock_vol_unit_out"])

        print("[DEBUG] worker cmd:", cmd)
        rc = shell.run(cmd, shell_mode=True)

        json_path = out_dir / f"{prefix}.json"
        if rc and json_path.is_file():
            return {"json": File(str(json_path))}
        if rc:
            raise RuntimeError(f"Dilution worker failed (rc={rc}).")
        if not json_path.is_file():
            raise FileNotFoundError(f"JSON not produced: {json_path}")

        return {"json": File(str(json_path))}
