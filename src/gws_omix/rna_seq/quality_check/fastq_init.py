#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fastqc_task.py  –  GWS task “FastQC”

• lit un tableau metadata (TSV) listant les FASTQ(.gz)
• saute TOUTES les lignes commençant par « # »
• exécute FastQC (multi-thread) sur tous les fichiers
• dépose les rapports HTML/ZIP dans   <working_dir>/fastqc_result
"""

import csv
from pathlib import Path
from typing import List

import pandas as pd

from gws_core import (
    ConfigParams,
    ConfigSpecs,
    Folder,          # classe générique de dossier
    InputSpec,
    InputSpecs,
    IntParam,
    OutputSpec,
    OutputSpecs,
    Task,
    TaskInputs,
    TaskOutputs,
    task_decorator,
    File,            # classe générique de fichier
)

from .fastq_init_env import FastqInitShellProxyHelper
from gws_omix import FastqFolder


# ───────────────────────────────  HELPERS  ────────────────────────────────
def _read_metadata(path: Path) -> pd.DataFrame:
    """
    Charger un TSV de métadonnées en ignorant toutes les lignes
    commençant par « # ». La première ligne non commentée devient l’en-tête.
    """
    return pd.read_csv(
        path,
        sep="\t",
        comment="#",          # ← saute entièrement les lignes #
        header=0,             # la 1ʳᵉ ligne restante = en-tête
        quoting=csv.QUOTE_NONE,
        dtype=str,
    )


def _normalise_paths(files: List[str], root: Path) -> List[str]:
    """Rendre tous les chemins absolus et vérifier qu’ils existent."""
    resolved: List[str] = []
    for p in files:
        path = Path(p)
        if not path.is_absolute():
            path = root / path
        if not path.exists():
            raise FileNotFoundError(f"Input file not found: {path}")
        resolved.append(str(path))
    return resolved


# ───────────────────────────────  TASK  ────────────────────────────────
@task_decorator("FastQC", human_name="FastQC",
                short_description="Check read quality with FastQC")
class FastqcInit(Task):

    # -------------  Spécifications des entrées / sorties  -------------
    input_specs: InputSpecs = InputSpecs(
        {
            "fastq_folder": InputSpec(
                FastqFolder,
                human_name="FASTQ folder",
                short_description="Folder containing *.fastq(.gz) files ; "
                                 "relative paths in metadata are resolved from here",
            ),
            "metadata": InputSpec(
                File,
                human_name="Metadata (TSV)",
                short_description="TSV listing sample-ids and file path(s)",
            ),
        }
    )

    output_specs: OutputSpecs = OutputSpecs(
        {
            "output": OutputSpec(
                Folder,
                human_name="Quality reports",
                short_description="Folder with FastQC HTML & ZIP files",
            )
        }
    )

    config_specs: ConfigSpecs = ConfigSpecs({
        "threads": IntParam(
            default_value=8,
            min_value=2,
            short_description="Number of CPU threads for FastQC",
        ),
    })

    # -------------------------------  RUN  -------------------------------
    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """Execute FastQC on all FASTQ files declared in *metadata*."""

        # -------- inputs & parameters
        fastq_root: FastqFolder = inputs["fastq_folder"]
        metadata_file: File = inputs["metadata"]
        threads: int = params["threads"]

        # -------- load metadata (skip every line starting with '#')
        df = _read_metadata(Path(metadata_file.path))

        paired_cols = {"forward-absolute-filepath", "reverse-absolute-filepath"}
        single_cols = {"absolute-filepath"}

        if paired_cols.issubset(df.columns):
            fastq_paths = (
                df["forward-absolute-filepath"].tolist()
                + df["reverse-absolute-filepath"].tolist()
            )
        elif single_cols.issubset(df.columns):
            fastq_paths = df["absolute-filepath"].tolist()
        else:
            raise ValueError(
                "Metadata must contain either 'absolute-filepath' (single-end) "
                "or 'forward-absolute-filepath' & 'reverse-absolute-filepath' (paired-end)."
            )

        # -------- make paths absolute & check existence
        fastq_paths = _normalise_paths(fastq_paths, Path(fastq_root.path))
        if not fastq_paths:
            raise ValueError("Metadata table contains no FASTQ paths.")

        # -------- working directory
        shell = FastqInitShellProxyHelper.create_proxy(self.message_dispatcher)
        result_dir = Path(shell.working_dir) / "fastqc_result"
        result_dir.mkdir(parents=True, exist_ok=True)

        # -------- run FastQC
        cmd = f"fastqc --noextract -t {threads} {' '.join(fastq_paths)} -o {result_dir}"
        if shell.run(cmd, shell_mode=True) != 0:
            raise RuntimeError("FastQC exited with non-zero status (see logs).")

        # -------- outputs
        return {"output": Folder(str(result_dir))}
