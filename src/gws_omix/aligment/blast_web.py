#!/usr/bin/env python3
# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
from pathlib import Path
from typing import Final

from gws_core import (
    ConfigParams, ConfigSpecs, File,
    FloatParam, InputSpec, InputSpecs,
    OutputSpec, OutputSpecs, ResourceSet,
    ShellProxy, StrParam, TableImporter,
    Task, TaskInputs, TaskOutputs, task_decorator
)

from ..base_env.blast_web_env_task import BlastWebShellProxyHelper


@task_decorator(
    "BlastWebRunner",
    human_name="NCBI_Remote_Blast",
    short_description="Run remote BLAST per sequence using NCBI's servers")

class BlastWebRunner(Task):

    input_specs: Final[InputSpecs] = InputSpecs({
        "input_fasta": InputSpec(File, human_name="FASTA file")
    })

    output_specs: Final[OutputSpecs] = OutputSpecs({
        "blast_results": OutputSpec(ResourceSet, human_name="BLAST Results")
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "blast_program": StrParam(
            allowed_values=["blastn", "blastp", "blastx", "tblastn"],
            short_description="BLAST program"
        ),
        "sequence_type": StrParam(
            allowed_values=["nucl", "prot"],
            short_description="Type of sequences in input"
        ),
        "evalue": FloatParam(default_value=1e-3, min_value=0.0),
        "max_target_seqs": FloatParam(default_value=20, min_value=1)
    })

    python_file_path: Final[str] = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_blast_web.py"
    )

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        fasta: File = inputs["input_fasta"]
        prog = params["blast_program"]
        seq_type = params["sequence_type"]
        evalue = params["evalue"]
        max_target_seqs = params["max_target_seqs"]
        pause = 2
        threads = 1

        shell: ShellProxy = BlastWebShellProxyHelper.create_proxy(self.message_dispatcher)

        output_dir = os.path.join(shell.working_dir, "blast_results")
        split_dir = os.path.join(shell.working_dir, "split_fasta")

        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(split_dir, exist_ok=True)

        cmd = (
            f"python3 {self.python_file_path} "
            f"{fasta.path} {output_dir} {split_dir} {seq_type} {prog} "
            f"{evalue} {int(max_target_seqs)} 1 {threads} {pause}"
        )

        if shell.run(cmd, shell_mode=True) != 0:
            raise RuntimeError("Remote BLAST failed")

        result_set = ResourceSet()
        for file in sorted(Path(output_dir).glob("*.tsv")):
            result_set.add_resource(
                TableImporter.call(File(str(file)), {
                    "delimiter": "tab",
                    "header": 0,
                    "file_format": "tsv"
                })
            )

        return {"blast_results": result_set}
