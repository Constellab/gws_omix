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
    "temperature_conversion",
    human_name="Temperature Conversion",
    short_description="Convert a temperature value between Celsius (°C), Fahrenheit (°F), and Kelvin (°K)."
)
class TemperatureConversion(Task):
    """
    Temperature Conversion
    -----------------------

    Converts a temperature value between:
      - **°C (Celsius)**
      - **°F (Fahrenheit)**
      - **°K (Kelvin)**

    **Formulas Used**
      - °C = (°F − 32) × 5/9
      - °F = (°C × 9/5) + 32
      - °K = °C + 273.16

    **Examples**
      - 25°C → 77°F → 298.16°K
      - 0°F → -17.78°C → 255.38°K
      - 310°K → 36.84°C → 98.33°F
    """

    input_specs: Final[InputSpecs] = InputSpecs({})
    output_specs: Final[OutputSpecs] = OutputSpecs({
        "json": OutputSpec(File, human_name="Result JSON")
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(default_value="temperature_conversion", short_description="Output JSON prefix"),
        "value": StrParam(default_value="", short_description="Temperature value to convert"),
        "unit": StrParam(
            default_value="C",
            allowed_values=["C", "F", "K"],
            short_description="Temperature unit: C (Celsius), F (Fahrenheit), K (Kelvin)"
        )
    })

    def _q(self, x) -> str:
        return shlex.quote(str(x))

    def _opt(self, name: str, val: str) -> str:
        v = (val or "").strip()
        return f" --{name} {self._q(v)}" if v else ""

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        shell: ShellProxy = DsdnaUgToPmolShellProxyHelper.create_proxy(self.message_dispatcher)
        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "temperature_conversion"
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix = (params["prefix"] or "temperature_conversion").strip() or "temperature_conversion"
        worker = Path(__file__).with_name("_temperature_conversion.py")
        if not worker.exists():
            raise FileNotFoundError(f"Worker not found: {worker}")

        cmd = f"python3 {self._q(worker)} --out {self._q(out_dir)} --prefix {self._q(prefix)}"
        cmd += self._opt("value", params["value"])
        cmd += self._opt("unit", params["unit"])

        print("[DEBUG] worker cmd:", cmd)
        rc = shell.run(cmd, shell_mode=True)

        json_path = out_dir / f"{prefix}.json"
        if rc and json_path.is_file():
            return {"json": File(str(json_path))}
        if rc:
            raise RuntimeError(f"Temperature Conversion worker failed (rc={rc}).")
        if not json_path.is_file():
            raise FileNotFoundError(f"JSON not produced: {json_path}")

        return {"json": File(str(json_path))}
