#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from typing import Final

from gws_core import (
    ConfigParams, ConfigSpecs, File,
    InputSpec, InputSpecs, OutputSpec, OutputSpecs,
    StrParam, Task, TaskInputs, TaskOutputs,
    task_decorator, ShellProxy,
)

from .dsdna_ug_to_pmol_env import DsdnaUgToPmolShellProxyHelper

@task_decorator(
    "dsDNA_ug_to_pmol",
    human_name="dsDNA: µg to pmol",
    short_description="Convertit une masse de dsDNA (µg) en pmol pour une longueur donnée (bp)."
)
class DSDNAUgToPmol(Task):
    """
    dsDNA: µg → pmol

    This module converts a mass of double-stranded DNA (µg) into a molecular quantity (pmol).
    In other words: how many molecules of a DNA fragment of a given length are contained in a given mass?

    Formula:
        pmol = (µg × 1e6) / (bp × 660)

    Where:
        - 1e6 = conversion factor from micrograms (µg) to picograms (pg)
        - 660 = average molecular weight of one base pair of double-stranded DNA (pg/pmol)
        - bp = length of the DNA fragment in base pairs

    Example:
        For 1 µg of a 3000 bp fragment:
            pmol = (1 × 1e6) / (3000 × 660) = 0.505 pmol

    Interpretation:
        → 0.505 pmol corresponds to approximately 3.04 × 10^11 DNA molecules
        (since 1 pmol = 6.022 × 10^11 molecules (number of Avogadro)).

    """
    input_specs: Final[InputSpecs] = InputSpecs({})
    output_specs: Final[OutputSpecs] = OutputSpecs({
        # 🔹 un seul fichier de sortie : le JSON
        "json": OutputSpec(File, human_name="Résultat JSON")
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(default_value="dsdna_ug_to_pmol", short_description="Préfixe du fichier JSON"),
        "length_bp": StrParam(default_value="", short_description="DNA Length (bp)"),
        "mass_ug": StrParam(default_value="", short_description="DNA Amount (µg)")
    })

    def _opt(self, name: str, val: str) -> str:
        v = (val or "").strip()
        return f' --{name} "{v}"' if v else ""

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        # Shell proxy fourni par votre infra (même que pour MAFFT)
        shell: ShellProxy = DsdnaUgToPmolShellProxyHelper.create_proxy(self.message_dispatcher)

        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "dsdna_ug_to_pmol"
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix = (params["prefix"] or "dsdna_ug_to_pmol").strip() or "dsdna_ug_to_pmol"
        worker = Path(__file__).with_name("_dsdna_ug_to_pmol.py")
        if not worker.exists():
            raise FileNotFoundError(f"Worker not found: {worker}")

        cmd = f'python3 "{worker}" --out "{out_dir}" --prefix "{prefix}"'
        cmd += self._opt("length-bp", params["length_bp"])
        cmd += self._opt("mass-ug", params["mass_ug"])

        print("[DEBUG] worker cmd:", cmd)
        rc = shell.run(cmd, shell_mode=True)

        json_path = out_dir / f"{prefix}.json"
        if rc:
            if json_path.is_file():
                return {"json": File(str(json_path))}
            raise RuntimeError(f"dsDNA_ug_to_pmol worker failed (rc={rc}).")

        if not json_path.is_file():
            raise FileNotFoundError(f"JSON not produced: {json_path}")

        return {"json": File(str(json_path))}
