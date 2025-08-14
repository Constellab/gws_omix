#!/usr/bin/env python3
# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import shutil
from pathlib import Path
from typing import Final

from gws_core import (
    ConfigParams, ConfigSpecs,
    File, Folder, InputSpec, InputSpecs,
    OutputSpec, OutputSpecs, ShellProxy,
    IntParam, FloatParam, StrParam,
    ResourceSet, TableImporter, Table,
    Task, TaskInputs, TaskOutputs,
    task_decorator
)

from gws_omix.base_env.blast_web_env_task import BlastWebShellProxyHelper


@task_decorator(
    "LocalBlast",
    human_name="RefSeqDB Local BLAST",
    short_description="Run local BLAST on each sequence of a FASTA using RefSeq DB")

class LocalBlast(Task):

    input_specs: Final[InputSpecs] = InputSpecs({
        "input_fasta": InputSpec(File, human_name="FASTA input"),
        "blast_db_folder": InputSpec(Folder, human_name="BLAST DB folder")
    })

    output_specs: Final[OutputSpecs] = OutputSpecs({
        "blast_results": OutputSpec(ResourceSet, human_name="BLAST output tables")
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "sequence_type": StrParam(allowed_values=["nucl", "prot"]),
        "blast_program": StrParam(allowed_values=["blastn", "blastp", "blastx", "tblastn"]),
        "evalue": FloatParam(default_value=1e-5),
        "max_target_seqs": IntParam(default_value=5),
        "threads": IntParam(default_value=8)
    })

    python_file_path: Final[str] = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_blast_local_task.py"
    )

    shell_proxy: ShellProxy = None
    split_dir: str = ""

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        fasta: File = inputs["input_fasta"]
        db_folder: Folder = inputs["blast_db_folder"]

        shell: ShellProxy = BlastWebShellProxyHelper.create_proxy(self.message_dispatcher)

        output_dir = os.path.join(shell.working_dir, "blast_results")
        self.split_dir = os.path.join(shell.working_dir, "split_fasta")
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(self.split_dir, exist_ok=True)

        cmd = (
            f"python3 {self.python_file_path} "
            f"{fasta.path} {params['sequence_type']} {params['blast_program']} "
            f"{db_folder.path} {output_dir} {self.split_dir} "
            f"{params['evalue']} {params['max_target_seqs']} 1 {params['threads']}"
        )

        print(f"[INFO] Running command: {cmd}")
        result = shell.run(cmd, shell_mode=True)
        if result != 0:
            raise RuntimeError("Local BLAST failed")

        self.shell_proxy = shell

        # collect result tables into a ResourceSet
        result_set = ResourceSet()
        for fn in sorted(os.listdir(output_dir)):
            if fn.endswith(".tsv"):
                path = os.path.join(output_dir, fn)
                table = TableImporter.call(
                    File(path),
                    {"delimiter": "tab", "header": 0, "file_format": "tsv"}
                )
                result_set.add_resource(table,fn)

        return {
            "blast_results": result_set
        }

    def run_after_task(self):
        if self.split_dir and os.path.exists(self.split_dir):
            shutil.rmtree(self.split_dir)

        if self.shell_proxy:
            self.shell_proxy.clean_working_dir()
