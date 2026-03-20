#!/usr/bin/env python3
# LICENSE
# Exclusive property of Gencovery SAS
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# Written by Gencovery SAS, 2024

from pathlib import Path
from typing import Final, List, Union

from gws_core import (
    ConfigParams,
    ConfigSpecs,
    File,
    Folder,
    InputSpec,
    InputSpecs,
    OutputSpec,
    OutputSpecs,
    ResourceSet,
    ShellProxy,
    Task,
    TaskInputs,
    TaskOutputs,
    task_decorator,
)

from .pycirclize_env import pyCirclizeShellProxyHelper


@task_decorator(
    "PyCirclize",
    human_name="Circular Genome Visualization (pyCirclize)",
    short_description="Circular visualization of GenBank or GFF3 genomes",
)
class PyCirclize(Task):
    """
    Visualizes one or more genome annotation files using pyCirclize.

    Input:
    - genome_input (File or Folder): GenBank/GFF3 file or folder

    Output:
    - circular_plots (ResourceSet): one PNG File per input genome file
    Community documentation : https://constellab.community/bricks/gws_omix/latest/doc/use-cases/circular-genome-visualization/
    """

    input_specs: Final[InputSpecs] = InputSpecs(
        {
            "genome_input": InputSpec(
                (Folder, File),
                human_name="Genome input",
                short_description=(
                    "A GenBank/GFF3 file or a folder containing one or more "
                    "supported files"
                ),
            )
        }
    )

    output_specs: Final[OutputSpecs] = OutputSpecs(
        {
            "circular_plots": OutputSpec(
                ResourceSet,
                human_name="Circular plots",
                short_description="One PNG per input genome file",
            ),
        }
    )

    config_specs: Final[ConfigSpecs] = ConfigSpecs({})

    def _list_dir_files(self, directory: Path) -> str:
        if not directory.is_dir():
            return "<directory does not exist>"

        names: List[str] = sorted(p.name for p in directory.iterdir())
        if not names:
            return "<no files>"

        return "\n".join(names)

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        genome_input: Union[Folder, File] = inputs["genome_input"]

        shell: ShellProxy = pyCirclizeShellProxyHelper.create_proxy(
            self.message_dispatcher
        )

        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "pycirclize"
        out_dir.mkdir(parents=True, exist_ok=True)

        worker = Path(__file__).with_name("_pyCirclize.py")
        input_path = genome_input.path

        cmd_parts = [
            "python3",
            str(worker),
            "--in",
            f"\"{input_path}\"",
            "--out",
            f"\"{out_dir}\"",
        ]
        cmd_str = " ".join(cmd_parts)

        rc = shell.run(cmd_str, shell_mode=True)

        if rc != 0:
            files_list = self._list_dir_files(out_dir)
            raise RuntimeError(
                f"pyCirclize worker failed with exit code {rc}.\n"
                f"Files present in {out_dir}:\n{files_list}"
            )

        png_files = sorted(out_dir.glob("*.png"))
        if not png_files:
            files_list = self._list_dir_files(out_dir)
            raise RuntimeError(
                f"pyCirclize completed but no PNG files were generated.\n"
                f"Files present in {out_dir}:\n{files_list}"
            )

        rs: ResourceSet = ResourceSet()
        rs.name = "Circular genome plots"
        for png in png_files:
            rs.add_resource(File(str(png)), png.name)

        return {
            "circular_plots": rs,
        }