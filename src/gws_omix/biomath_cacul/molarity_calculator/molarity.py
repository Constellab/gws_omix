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
    "molarity_calculator",
    human_name="Molarity Calculator",
    short_description="Compute mass to weigh from MW, target molarity, and final volume."
)
class MolarityCalculator(Task):
    """
    Molarity Calculator
    -------------------

    Calculates the mass of a compound to weigh given:
      - Molecular Weight (g/mol)
      - Desired final concentration (molarity)
      - Desired final volume

    Relations:
      moles = M × V
      mass (g) = moles × MW = M × V × MW

    """

    input_specs: Final[InputSpecs] = InputSpecs({})
    output_specs: Final[OutputSpecs] = OutputSpecs({
        "json": OutputSpec(File, human_name="Result JSON")
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(default_value="molarity", short_description="Output JSON prefix"),

        "mw": StrParam(default_value="", short_description="Molecular Weight (g/mol)"),

        "molarity": StrParam(default_value="", short_description="Final concentration"),
        "mol_unit": StrParam(
            default_value="mM",
            allowed_values=["M", "mM", "µM"],
            short_description="Unit for final concentration"
        ),

        "final_vol": StrParam(default_value="", short_description="Final volume"),
        "vol_unit": StrParam(
            default_value="mL",
            allowed_values=["L", "mL", "µL"],
            short_description="Unit for final volume"
        ),

        "mass_unit_out": StrParam(
            default_value="g",
            allowed_values=["g", "mg", "µg"],
            short_description="Output mass unit."
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
        out_dir = work_dir / "molarity"
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix = (params["prefix"] or "molarity").strip() or "molarity"
        worker = Path(__file__).with_name("_molarity.py")
        if not worker.exists():
            raise FileNotFoundError(f"Worker not found: {worker}")

        cmd = f"python3 {self._q(worker)} --out {self._q(out_dir)} --prefix {self._q(prefix)}"
        cmd += self._opt("mw", params["mw"])
        cmd += self._opt("molarity", params["molarity"])
        cmd += self._opt("mol-unit", params["mol_unit"])
        cmd += self._opt("final-vol", params["final_vol"])
        cmd += self._opt("vol-unit", params["vol_unit"])
        cmd += self._opt("mass-unit-out", params["mass_unit_out"])

        print("[DEBUG] worker cmd:", cmd)
        rc = shell.run(cmd, shell_mode=True)

        json_path = out_dir / f"{prefix}.json"
        if rc and json_path.is_file():
            return {"json": File(str(json_path))}
        if rc:
            raise RuntimeError(f"Molarity worker failed (rc={rc}).")
        if not json_path.is_file():
            raise FileNotFoundError(f"JSON not produced: {json_path}")

        return {"json": File(str(json_path))}
