#!/usr/bin/env python3
# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
from typing import Final

from gws_core import (
    ConfigParams,
    Folder,
    OutputSpec,
    OutputSpecs,
    ShellProxy,
    Task,
    TaskInputs,
    TaskOutputs,
    task_decorator,
)

from gws_omix.base_env.blast_web_env_task import BlastWebShellProxyHelper


@task_decorator(
    "Build_RefSeqBlastDB",
    human_name="Build_RefSeqBlast_Database",
    short_description="Download and prepare RefSeq BLAST databases (RNA & Protein)")

class RefSeqBlastDB(Task):

    input_specs = {}
    output_specs: Final[OutputSpecs] = OutputSpecs({
        "blast_db_folder": OutputSpec(Folder, human_name="RefSeq BLAST DB folder")
    })

    python_file_path: Final[str] = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_blast_refseq_makedb.py"
    )

    shell_proxy: ShellProxy = None

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        shell: ShellProxy = BlastWebShellProxyHelper.create_proxy(self.message_dispatcher)

        output_dir = os.path.join(shell.working_dir, "refseq_blastdb")
        os.makedirs(output_dir, exist_ok=True)

        cmd = f"python3 {self.python_file_path} {output_dir}"
        print(f"[INFO] Running shell command: {cmd}")

        result = shell.run(cmd, shell_mode=True)
        if result != 0:
            raise RuntimeError("RefSeq BLAST DB creation failed")

        self.shell_proxy = shell
        return {
            "blast_db_folder": Folder(output_dir)
        }
