# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import csv, os
from pathlib import Path
from typing import List

import pandas as pd
from gws_core import ( ConfigParams, ConfigSpecs, Folder, File, InputSpec, InputSpecs, OutputSpec, OutputSpecs,
                        IntParam, Task, TaskInputs, TaskOutputs, task_decorator)
                        
from .hisat2_env import Hisat2ShellProxyHelper
from gws_omix import FastqFolder


def _read_meta(fp: Path) -> pd.DataFrame:
    return pd.read_csv(fp, sep="\t", comment="#", header=0,
                       quoting=csv.QUOTE_NONE, dtype=str)


def _abs(p: str, root: Path) -> str:
    q = Path(p) if Path(p).is_absolute() else root / p
    if not q.exists():
        raise FileNotFoundError(q)
    return str(q)


@task_decorator("Hisat2_Align", human_name="Hisat2_Align",
                short_description="Align RNA-seq reads with HISAT2 (keeps Sample names)")
class Hisat2Align(Task):

    input_specs: InputSpecs = InputSpecs(
        {
            "fastq_folder": InputSpec(
                FastqFolder, human_name="FASTQ folder",
                short_description="Base folder for relative paths"),
            "metadata": InputSpec(
                File, human_name="Metadata (TSV)",
                short_description="Must contain Sample column + file path columns"),
            "hisat2_index": InputSpec(
                Folder, human_name="HISAT2 index",
                short_description="Folder containing genome_index.*.ht2"),
        }
    )
    output_specs: OutputSpecs = OutputSpecs(
        {"output": OutputSpec(Folder, human_name="BAM folder")}
    )
    config_specs: ConfigSpecs = ConfigSpecs({
        "threads": IntParam(default_value=4, min_value=2),
    })

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:

        fastq_root: FastqFolder = inputs["fastq_folder"]
        meta_file : File   = inputs["metadata"]
        index_dir : Folder = inputs["hisat2_index"]

        root = Path(fastq_root.path)
        index_prefix = Path(index_dir.path) / "genome_index"
        threads = params["threads"]

        meta = _read_meta(Path(meta_file.path))

        is_pe = {"forward-absolute-filepath", "reverse-absolute-filepath"} <= set(meta.columns)
        is_se = "absolute-filepath" in meta.columns
        if not (is_pe or is_se):
            raise ValueError("Metadata lacks required filepath columns.")

        shell = Hisat2ShellProxyHelper.create_proxy(self.message_dispatcher)
        out_dir = Path(shell.working_dir) / "result"
        out_dir.mkdir(parents=True, exist_ok=True)

        if is_se:
            print("[INFO] single-end alignment")
            for _, row in meta.iterrows():
                sample = row.get("Sample", row.get("sample-id", "sample")).strip()
                fq = _abs(row["absolute-filepath"], root)
                bam = out_dir / f"{sample}.bam"

                cmd = (
                    f"hisat2 -q -p {threads} "
                    f"-x {index_prefix} -U {fq} | "
                    f"samtools sort -@ {threads} -o {bam}"
                )
                print("[DEBUG]", cmd)
                if shell.run(cmd, shell_mode=True):
                    raise RuntimeError(f"HISAT2 failed on {sample}")

        else:
            print("[INFO] paired-end alignment")
            for _, row in meta.iterrows():
                sample = row.get("Sample", row.get("sample-id", "sample")).strip()
                fq1 = _abs(row["forward-absolute-filepath"],  root)
                fq2 = _abs(row["reverse-absolute-filepath"], root)
                bam = out_dir / f"{sample}.bam"

                cmd = (
                    f"hisat2 -q -p {threads} "
                    f"-x {index_prefix} -1 {fq1} -2 {fq2} | "
                    f"samtools sort -@ {threads} -o {bam}"
                )
                print("[DEBUG]", cmd)
                if shell.run(cmd, shell_mode=True):
                    raise RuntimeError(f"HISAT2 failed on {sample}")

        return {"output": Folder(str(out_dir))}
