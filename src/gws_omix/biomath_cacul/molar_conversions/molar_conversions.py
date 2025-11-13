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
    "molar_conversions",
    human_name="Molar Conversions (Protein)",
    short_description="Convert between Protein Size (kDa), pmol of Protein, and µg of Protein."
)
class MolarConversions(Task):
    """
    Molar Conversions

    This module converts between:
    - Protein Size (kDa)
    - pmol of Protein
    - µg of Protein

    - µg protein = protein size (kDa) × pmol protein × (10^9 µg / kg) × (kg / 10^12 pmol)
    - pmol protein = (µg protein / protein size (kDa)) × (10^12 pmol / mol) × (mol / 10^9 µg)
    - protein size (kDa) = (µg protein / pmol protein) × (10^12 pmol / mol) × (mol / 10^9 µg)

    Note: A protein that is 1 kDa has a molecular weight of 1 kilogram per mole.

    Derivations (step-by-step)
    ==========================

    1) µg of Protein
    ----------------
    MW in µg/mol = kDa × (1 kg/mol) × (10^9 µg / 1 kg) = kDa × 10^9 µg/mol
    moles        = pmol × (10^-12 mol / 1 pmol)

    µg = (MW in µg/mol) × (moles)
    = (kDa × 10^9 µg/mol) × (pmol × 10^-12 mol)
    = kDa × pmol × 10^-3
    = (kDa × pmol) / 1000

    This matches the Promega line that uses “×10^9 µg/kg × mol / 10^12 pmol” (i.e., kg→µg and mol→pmol conversions).


    2) pmol of Protein
    ------------------
    pmol = (µg / (MW in µg/mol)) × 10^12
        = (µg / (kDa × 10^9)) × 10^12
        = (µg × 10^3) / kDa

    So: pmol = (µg × 1000) / kDa


    3) Protein Size (kDa)
    ---------------------
    kDa = MW / (10^9 µg/mol)
    But MW = (µg / pmol) × 10^12 (since pmol × 10^-12 = mol → µg per mol is scaled accordingly)

    Therefore:
    kDa = [ (µg / pmol) × 10^12 ] / 10^9
        = (µg × 10^3) / pmol

    So: kDa = (µg × 1000) / pmol


    Simplified formulas used by this tool
    =====================================
    - µg  = (kDa × pmol) / 1000
    - pmol = (µg × 1000) / kDa
    - kDa  = (µg × 1000) / pmol

        If you choose: Protein Size (kDa)
      • Fill in: µg of Protein, pmol of Protein
      • Formula: kDa = (µg × 1000) / pmol

    If you choose: µg of Protein
      • Fill in: Protein Size (kDa), pmol of Protein
      • Formula: µg = (kDa × pmol) / 1000

    If you choose: pmol of Protein
      • Fill in: Protein Size (kDa), µg of Protein
      • Formula: pmol = (µg × 1000) / kDa
    """


    input_specs: Final[InputSpecs] = InputSpecs({})
    output_specs: Final[OutputSpecs] = OutputSpecs({
        "json": OutputSpec(File, human_name="Result JSON")
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(default_value="molar_conversions", short_description="Output JSON prefix"),
        "calculate_for": StrParam(
            default_value="Protein Size (kDa)",
            allowed_values=["Protein Size (kDa)", "pmol of Protein", "μg of Protein"],
            short_description=(
                "Select what to calculate. Requirements:\n"
                "• Protein Size (kDa): requires µg of Protein AND pmol of Protein\n"
                "• pmol of Protein: requires Protein Size (kDa) AND µg of Protein\n"
                "• µg of Protein: requires Protein Size (kDa) AND pmol of Protein"
            ),
        ),
        "size_kda": StrParam(default_value="", short_description="Protein Size (kDa)"),
        "pmol": StrParam(default_value="", short_description="pmol of Protein"),
        "mass_ug": StrParam(default_value="", short_description="µg of Protein"),
    })

    def _q(self, x) -> str:
        """Shell-safe quote for one argument."""
        return shlex.quote(str(x))

    def _opt(self, name: str, val: str) -> str:
        v = (val or "").strip()
        return f" --{name} {self._q(v)}" if v else ""

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        shell: ShellProxy = DsdnaUgToPmolShellProxyHelper.create_proxy(self.message_dispatcher)

        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "molar_conversions"
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix = (params["prefix"] or "molar_conversions").strip() or "molar_conversions"
        worker = Path(__file__).with_name("_molar_conversions.py")
        if not worker.exists():
            raise FileNotFoundError(f"Worker not found: {worker}")

        # Build a shell-safe command without inner quotes
        cmd = f"python3 {self._q(worker)} --out {self._q(out_dir)} --prefix {self._q(prefix)}"
        cmd += self._opt("calculate-for", params["calculate_for"])
        cmd += self._opt("size-kda", params["size_kda"])
        cmd += self._opt("pmol", params["pmol"])
        cmd += self._opt("mass-ug", params["mass_ug"])

        print("[DEBUG] worker cmd:", cmd)
        rc = shell.run(cmd, shell_mode=True)

        json_path = out_dir / f"{prefix}.json"
        if rc and json_path.is_file():
            return {"json": File(str(json_path))}
        if rc:
            raise RuntimeError(f"Molar Conversions worker failed (rc={rc}).")
        if not json_path.is_file():
            raise FileNotFoundError(f"JSON not produced: {json_path}")

        return {"json": File(str(json_path))}
