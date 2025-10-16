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
    "OD260_to_ugml",
    human_name="Nucleic Acid: OD260 to µg/ml",
    short_description="Convert OD260 to µg/ml for DNA/RNA/ssDNA/ssOligo"
)
class OD260ToUgml(Task):

    """
    Nucleic Acid: OD260 → µg/ml

    Formula:
        µg/ml = OD260 × conversion factor

    Conversion factors:
        - DNA (dsDNA): 50 µg/ml per OD
        - RNA (ssRNA): 40 µg/ml per OD
        - Single-stranded DNA: 35 µg/ml per OD
        - Single-stranded Oligo: 20 µg/ml per OD
    """

    input_specs: Final[InputSpecs] = InputSpecs({})
    output_specs: Final[OutputSpecs] = OutputSpecs({
        "json": OutputSpec(File, human_name="Résultat JSON")
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(default_value="od260_to_ugml", short_description="Préfixe du fichier JSON"),
        "od260": StrParam(default_value="", short_description="Valeur OD à 260 nm"),
        "sample_type": StrParam(
            default_value="DNA",
            allowed_values=["DNA", "RNA", "Single-stranded-DNA", "Single-stranded-Oligo"],
            short_description="Type d'échantillon : DNA | RNA | Single-stranded-DNA | Single-stranded-Oligo"
        )
    })

    def _opt(self, name: str, val: str) -> str:
        v = (val or "").strip()
        return f' --{name} "{v}"' if v else ""

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        shell: ShellProxy = DsdnaUgToPmolShellProxyHelper.create_proxy(self.message_dispatcher)
        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "od260_to_ugml"
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix = params["prefix"].strip() or "od260_to_ugml"
        worker = Path(__file__).with_name("_od260_to_ugml.py")
        if not worker.exists():
            raise FileNotFoundError(f"Worker introuvable : {worker}")

        cmd = f'python3 "{worker}" --out "{out_dir}" --prefix "{prefix}"'
        cmd += self._opt("od260", params["od260"])
        cmd += self._opt("sample-type", params["sample_type"])

        print("[DEBUG] worker cmd:", cmd)
        rc = shell.run(cmd, shell_mode=True)

        json_path = out_dir / f"{prefix}.json"
        if rc and json_path.is_file():
            return {"json": File(str(json_path))}
        if rc:
            raise RuntimeError(f"OD260_to_ugml worker failed (rc={rc}).")
        if not json_path.is_file():
            raise FileNotFoundError(f"JSON non produit : {json_path}")

        return {"json": File(str(json_path))}
