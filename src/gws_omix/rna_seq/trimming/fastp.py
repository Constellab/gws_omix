# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import csv
from pathlib import Path
from typing import List

import pandas as pd
from gws_core import (
    ConfigParams,
    ConfigSpecs,
    Folder, File,
    InputSpec, InputSpecs,
    OutputSpec, OutputSpecs,
    IntParam, StrParam,
    Task, TaskInputs, TaskOutputs,
    task_decorator,
)
from .fastp_env import FastpShellProxyHelper
from gws_omix import FastqFolder

def _read_metadata(path: Path) -> pd.DataFrame:
    """Read TSV; skip comment lines (#)."""
    return pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=0,
        quoting=csv.QUOTE_NONE,
        dtype=str,
    )


def _abs(p: str, root: Path) -> str:
    """Absolute path anchored on *root*, ensure it exists."""
    fp = Path(p) if Path(p).is_absolute() else root / p
    if not fp.exists():
        raise FileNotFoundError(fp)
    return str(fp)


@task_decorator(
    "Fastp",
    human_name="Fastp",
    short_description="Adapter & quality trimming with fastp (metadata driven)",
)
class Fastp(Task):

    input_specs: InputSpecs = InputSpecs(
        {
            "fastq_folder": InputSpec(
                FastqFolder,
                human_name="FASTQ folder",
                short_description="Folder containing raw FASTQ(.gz); "
                                 "relative paths in metadata start from here",
            ),
            "metadata": InputSpec(
                File,
                human_name="Metadata (TSV)",
                short_description="Table with absolute-filepath (SE) or "
                                 "forward- & reverse-absolute-filepath (PE)",
            ),
        }
    )

    output_specs: OutputSpecs = OutputSpecs(
        {
            "output": OutputSpec(
                FastqFolder,
                human_name="Trimmed FASTQ",
                short_description="Folder with *_trimmed*.fastq.gz files",
            )
        }
    )

    config_specs: ConfigSpecs = ConfigSpecs({
        "threads": IntParam(default_value=4, min_value=2),
        "Forward_separator": StrParam(default_value="1"),
        "Reverse_separator": StrParam(default_value="2"),
        "5_prime_hard_trimming_read_size": IntParam(default_value=0, min_value=0),
    })

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:

        fastq_root: FastqFolder = inputs["fastq_folder"]
        meta_file : File   = inputs["metadata"]

        root_path = Path(fastq_root.path)
        meta      = _read_metadata(Path(meta_file.path))

        threads  = params["threads"]
        fwd_sep  = params["Forward_separator"]
        rev_sep  = params["Reverse_separator"]
        headcrop = params["5_prime_hard_trimming_read_size"]

        is_paired = {"forward-absolute-filepath", "reverse-absolute-filepath"} <= set(meta.columns)
        is_single = "absolute-filepath" in meta.columns

        if not (is_single or is_paired):
            raise ValueError(
                "Metadata must contain either 'absolute-filepath' (single-end) "
                "or 'forward-absolute-filepath' & 'reverse-absolute-filepath' (paired-end)."
            )

        shell = FastpShellProxyHelper.create_proxy(self.message_dispatcher)
        result_dir = Path(shell.working_dir) / "result"
        result_dir.mkdir(parents=True, exist_ok=True)

        if is_single:
            print("[INFO] Single-end trimming with fastp")
            for _, row in meta.iterrows():
                sample = row.get("Sample", row.get("sample-id", "sample")).strip()
                in_fp  = _abs(row["absolute-filepath"], root_path)
                out_fp = result_dir / f"{sample}.fastq.gz"

                cmd = (
                    f"fastp --in1 {in_fp} --out1 {out_fp} "
                    f"--trim_front1 {headcrop} "
                    f"--thread {threads}"
                )
                print("[DEBUG]", cmd)
                if shell.run(cmd, shell_mode=True) != 0:
                    raise RuntimeError(f"fastp failed on sample {sample}")

        else:
            print("[INFO] Paired-end trimming with fastp")
            for _, row in meta.iterrows():
                sample = row.get("Sample", row.get("sample-id", "sample")).strip()
                fwd_in = _abs(row["forward-absolute-filepath"],  root_path)
                rev_in = _abs(row["reverse-absolute-filepath"], root_path)

                out1 = result_dir / f"{sample}_{fwd_sep}.fastq.gz"
                out2 = result_dir / f"{sample}_{rev_sep}.fastq.gz"

                cmd = (
                    f"fastp --in1 {fwd_in} --in2 {rev_in} "
                    f"--out1 {out1} --out2 {out2} "
                    f"--trim_front1 {headcrop} --trim_front2 {headcrop} "
                    f"--detect_adapter_for_pe "
                    f"--thread {threads}"
                )
                print("[DEBUG]", cmd)
                if shell.run(cmd, shell_mode=True) != 0:
                    raise RuntimeError(f"fastp failed on sample {sample}")

        return {"output": FastqFolder(str(result_dir))}
