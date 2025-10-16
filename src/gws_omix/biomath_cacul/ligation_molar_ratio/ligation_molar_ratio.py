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
from gws_omix.base_env.msaviz_env_task import MSAVisShellProxyHelper


@task_decorator(
    "Ligation_Molar_Ratio",
    human_name="Ligations: Molar Ratio of Insert:Vector",
    short_description="Compute insert mass (ng) for a chosen I:V molar ratio."
)
class LigationMolarRatio(Task):
    """
Ligations: Molar Ratio of Insert:Vector

This calculator determines how many nanograms (ng) of an **insert DNA fragment**
are required to achieve a specific **molar ratio** with a given vector DNA fragment.

Formula:
    ng_insert = (I / V) × (kb_insert / kb_vector) × ng_vector

Where:
    - I = desired molar ratio of insert (e.g., 1 for 1:1, or 3 for 3:1)
    - V = molar ratio of vector (usually 1)
    - kb_insert = insert length in kilobases (kb)
    - kb_vector = vector length in kilobases (kb)
    - ng_vector = amount of vector DNA used in nanograms (ng)

Example:
    For a 0.4 kb insert, 0.744 kb vector, 8 ng of vector DNA, and a 1:1 ratio:
        ng_insert = (1 / 1) × (0.4 / 0.744) × 8 = 4.3 ng
        
    """
    input_specs: Final[InputSpecs] = InputSpecs({})
    output_specs: Final[OutputSpecs] = OutputSpecs({
        "json": OutputSpec(File, human_name="Result JSON")
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(default_value="ligation_molar_ratio", short_description="Output JSON prefix"),
        "insert_kb": StrParam(default_value="", short_description="Insert Length (kb)"),
        "vector_ng": StrParam(default_value="", short_description="Vector Amount (ng)"),
        "vector_kb": StrParam(default_value="", short_description="Vector Length (kb)"),
        "ratio_insert": StrParam(default_value="1", short_description="Insert ratio I (default 1)"),
        "ratio_vector": StrParam(default_value="1", short_description="Vector ratio V (default 1)")
    })

    def _opt(self, name: str, val: str) -> str:
        v = (val or "").strip()
        return f' --{name} "{v}"' if v else ""

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        shell: ShellProxy = MSAVisShellProxyHelper.create_proxy(self.message_dispatcher)

        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "ligation_molar_ratio"
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix = (params["prefix"] or "ligation_molar_ratio").strip() or "ligation_molar_ratio"
        worker = Path(__file__).with_name("_ligation_molar_ratio.py")
        if not worker.exists():
            raise FileNotFoundError(f"Worker not found: {worker}")

        cmd = f'python3 "{worker}" --out "{out_dir}" --prefix "{prefix}"'
        cmd += self._opt("insert-kb", params["insert_kb"])
        cmd += self._opt("vector-ng", params["vector_ng"])
        cmd += self._opt("vector-kb", params["vector_kb"])
        cmd += self._opt("ratio-insert", params["ratio_insert"])
        cmd += self._opt("ratio-vector", params["ratio_vector"])

        print("[DEBUG] worker cmd:", cmd)
        rc = shell.run(cmd, shell_mode=True)

        json_path = out_dir / f"{prefix}.json"
        if rc and json_path.is_file():
            return {"json": File(str(json_path))}
        if rc:
            raise RuntimeError(f"Ligation_Molar_Ratio worker failed (rc={rc}).")

        if not json_path.is_file():
            raise FileNotFoundError(f"JSON not produced: {json_path}")

        return {"json": File(str(json_path))}
