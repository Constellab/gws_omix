#!/usr/bin/env python3
# LICENSE
# Exclusive property of Gencovery SAS
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# Written by Gencovery SAS, 2024

from pathlib import Path
from typing import Final, List

from gws_core import (
    ConfigParams, ConfigSpecs, File, Folder,
    InputSpec, InputSpecs, OutputSpec, OutputSpecs,
    ShellProxy, StrParam, Task, TaskInputs, TaskOutputs,
    task_decorator
)

from .pygenomeviz_env import PyGenomeVizShellProxyHelper


@task_decorator(
    "PyGenomeViz",
    human_name="Genome Visualization (pyGenomeViz)",
    short_description="Interactive visualization of multiple GenBank genomes"
)
class PyGenomeViz(Task):
    """
    Visualizes multiple GenBank genomes using pyGenomeViz core API.

    **Input**
    - `genbank_folder` (Folder): Folder containing one or more GenBank files (.gb, .gbk, .genbank)

    **Output**
    - `interactive_html` (File): The generated interactive HTML visualization
    """

    input_specs: Final[InputSpecs] = InputSpecs({
        "genbank_folder": InputSpec(
            Folder,
            human_name="GenBank folder",
            short_description="Folder containing one or more GenBank files"
        )
    })

    output_specs: Final[OutputSpecs] = OutputSpecs({
        "interactive_html": OutputSpec(
            File,
            human_name="Interactive HTML output"
        )
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(
            default_value="pygenomeviz",
            short_description="Output HTML file prefix"
        )
    })

    def _list_dir_files(self, directory: Path) -> str:
        if not directory.is_dir():
            return "<directory does not exist>"
        names: List[str] = sorted(p.name for p in directory.iterdir())
        if not names:
            return "<no files>"
        return "\n".join(names)

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        genbank_folder: Folder = inputs["genbank_folder"]

        shell: ShellProxy = PyGenomeVizShellProxyHelper.create_proxy(self.message_dispatcher)

        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "pygenomeviz"
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix = (params["prefix"] or "pygenomeviz").strip() or "pygenomeviz"
        worker = Path(__file__).with_name("_pyGenomeViz.py")

        expected_html = out_dir / f"{prefix}.html"

        cmd_parts = [
            "python3",
            str(worker),
            "--in", f"\"{genbank_folder.path}\"",
            "--out", f"\"{out_dir}\"",
            "--prefix", f"\"{prefix}\"",
        ]
        cmd_str = " ".join(cmd_parts)
        print("[DEBUG] pyGenomeViz worker cmd:", cmd_str)

        rc = shell.run(cmd_str, shell_mode=True)

        if rc != 0:
            files_list = self._list_dir_files(out_dir)
            raise RuntimeError(
                f"pyGenomeViz worker failed with exit code {rc}.\n"
                f"Files present in {out_dir}:\n{files_list}\n"
            )

        if not expected_html.is_file():
            files_list = self._list_dir_files(out_dir)
            raise RuntimeError(
                "pyGenomeViz completed but the expected HTML file was not found.\n"
                f"Expected: {expected_html}\n"
                f"Files present in {out_dir}:\n{files_list}\n"
            )

        return {
            "interactive_html": File(str(expected_html))
        }