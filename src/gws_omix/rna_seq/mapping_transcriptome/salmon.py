# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import csv
import os
from pathlib import Path
from typing import List

import pandas as pd
from gws_core import (ConfigParams,ConfigSpecs,Folder,File,InputSpec,InputSpecs,OutputSpec,OutputSpecs,
                       IntParam,Task,TaskInputs,TaskOutputs,task_decorator)

from .salmon_env import SalmonShellProxyHelper
from gws_omix import FastqFolder


def _read_metadata(path: Path) -> pd.DataFrame:
    """Read TSV and ignore every line starting with ‘#’."""
    return pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=0,
        quoting=csv.QUOTE_NONE,
        dtype=str,
    )


def _abs_paths(files: List[str], root: Path) -> List[str]:
    """Return absolute, verified paths (root-relative if needed)."""
    out: List[str] = []
    for p in files:
        fp = Path(p) if Path(p).is_absolute() else root / p
        if not fp.exists():
            raise FileNotFoundError(fp)
        out.append(str(fp))
    return out


@task_decorator(
    "Salmon_Quant",
    human_name="Salmon_Quant",
    short_description="Quantify RNA-seq reads with Salmon (metadata driven)",
)
class SalmonQuant(Task):
    """
    This task quantifies RNA-seq reads with Salmon using a metadata-driven workflow.
    Provide a fastq_folder (base path for resolving relative file paths), a TSV metadata containing a Sample column plus either single-end (absolute-filepath) or paired-end (forward-absolute-filepath and reverse-absolute-filepath) columns, and a salmon_index built beforehand.
    For each row, the task runs salmon quant (with --validateMappings) and writes outputs to working_dir (including quant.sf and logs).
    Configure threads to control CPU usage.
    """
    input_specs: InputSpecs = InputSpecs(
        {
            "fastq_folder": InputSpec(
                FastqFolder,
                human_name="FASTQ folder",
                short_description="Folder containing FASTQ(.gz) files; "
                                 "relative paths in metadata are resolved from here",
            ),
            "metadata": InputSpec(
                File,
                human_name="Metadata (TSV)",
                short_description="Table with absolute-filepath (SE) or "
                                 "forward- & reverse-absolute-filepath (PE)",
            ),
            "salmon_index": InputSpec(
                Folder,
                human_name="Salmon index",
                short_description="Transcriptome index created with `salmon index`",
            ),
        }
    )

    output_specs: OutputSpecs = OutputSpecs(
        {
            "output": OutputSpec(
                Folder,
                human_name="Salmon quant folder",
                short_description="One sub-folder per sample with quant.sf, logs …",
            )
        }
    )

    config_specs: ConfigSpecs = ConfigSpecs({
        "threads": IntParam(default_value=4, min_value=2, short_description="CPU threads"),
    })

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:

        fastq_root: FastqFolder = inputs["fastq_folder"]
        meta_file: File = inputs["metadata"]
        salmon_index: Folder = inputs["salmon_index"]

        root_path = Path(fastq_root.path)
        index_path = salmon_index.path
        threads = params["threads"]

        meta = _read_metadata(Path(meta_file.path))
        is_paired = {"forward-absolute-filepath", "reverse-absolute-filepath"} <= set(meta.columns)
        is_single = "absolute-filepath" in meta.columns

        if not (is_single or is_paired):
            raise ValueError(
                "Metadata must contain :\n"
                "  • absolute-filepath (single-end)  OR\n"
                "  • forward-absolute-filepath & reverse-absolute-filepath (paired-end)"
            )

        if "Sample" not in meta.columns:
            raise ValueError("Metadata file must include a 'Sample' column for naming output directories.")

        shell = SalmonShellProxyHelper.create_proxy(self.message_dispatcher)
        result_dir = Path(shell.working_dir) / "result"
        result_dir.mkdir(parents=True, exist_ok=True)

        if is_single:
            print("[INFO] Detected single-end metadata")
            fastq_files = _abs_paths(meta["absolute-filepath"].tolist(), root_path)

            for sample, in_fp in zip(meta["Sample"], fastq_files):
                out_dir = result_dir / sample
                out_dir.mkdir(exist_ok=True)

                cmd = (
                    f"salmon quant -i {index_path} -p {threads} -l A "
                    f"-r {in_fp} --validateMappings -o {out_dir}"
                )
                print("[DEBUG]", cmd)
                if shell.run(cmd, shell_mode=True) != 0:
                    raise RuntimeError(f"Salmon failed (SE) on {sample}")

        else:
            print("[INFO] Detected paired-end metadata")
            fwd_files = _abs_paths(meta["forward-absolute-filepath"].tolist(), root_path)
            rev_files = _abs_paths(meta["reverse-absolute-filepath"].tolist(), root_path)

            for sample, fwd_fp, rev_fp in zip(meta["Sample"], fwd_files, rev_files):
                out_dir = result_dir / sample
                out_dir.mkdir(exist_ok=True)

                cmd = (
                    f"salmon quant -i {index_path} -p {threads} -l A "
                    f"-1 {fwd_fp} -2 {rev_fp} "
                    f"--validateMappings -o {out_dir}"
                )
                print("[DEBUG]", cmd)
                if shell.run(cmd, shell_mode=True) != 0:
                    raise RuntimeError(f"Salmon failed (PE) on {sample}")

        return {"output": Folder(str(result_dir))}
