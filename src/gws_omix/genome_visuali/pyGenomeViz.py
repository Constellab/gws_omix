#!/usr/bin/env python3
# LICENSE
# Exclusive property of Gencovery SAS
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# Written by Gencovery SAS, 2024

from pathlib import Path
from typing import Final, List

from gws_core import (
    ConfigParams,
    ConfigSpecs,
    File,
    Folder,
    InputSpec,
    InputSpecs,
    OutputSpec,
    OutputSpecs,
    ShellProxy,
    StrParam,
    Task,
    TaskInputs,
    TaskOutputs,
    task_decorator,
)

from .pygenomeviz_env import PyGenomeVizShellProxyHelper

COMPARISON_MODES = [
    "visualization_only",
    "genome_comparison",
    "cds_protein_homology",
]


@task_decorator(
    "PyGenomeVizCompare",
    human_name="Comparative Genome Visualization",
    short_description="Comparative genome visualization using pyGenomeViz",
)
class PyGenomeVizCompare(Task):
    
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
            "report_html": OutputSpec(
                File,
                human_name="HTML report",
            )
        }
    )

    config_specs: Final[ConfigSpecs] = ConfigSpecs(
        {
            "prefix": StrParam(
                default_value="pygenomeviz",
                short_description="Output file prefix",
            ),
            "comparison_mode": StrParam(
                allowed_values=COMPARISON_MODES,
                default_value="visualization_only",
                short_description=(
                    "visualization_only: display annotations only; "
                    "genome_comparison: compare whole genomes; "
                    "cds_protein_homology: compare homologous CDS/proteins."
                ),
            ),
        }
    )

    def _list_dir_files(self, directory: Path) -> str:
        if not directory.is_dir():
            return "<directory does not exist>"
        names: List[str] = sorted(p.name for p in directory.iterdir())
        return "\n".join(names) if names else "<no files>"

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        genome_input = inputs["genome_input"]
        shell: ShellProxy = PyGenomeVizShellProxyHelper.create_proxy(self.message_dispatcher)

        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "pygenomeviz"
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix = (params["prefix"] or "pygenomeviz").strip() or "pygenomeviz"
        comparison_mode = (params["comparison_mode"] or "visualization_only").strip()

        if comparison_mode not in COMPARISON_MODES:
            raise RuntimeError(
                f"Invalid comparison_mode: {comparison_mode}. "
                f"Expected one of: {', '.join(COMPARISON_MODES)}"
            )

        input_path = getattr(genome_input, "path", None)
        if not input_path:
            raise RuntimeError("Input does not provide a valid filesystem path.")

        worker = Path(__file__).with_name("_pyGenomeViz.py")
        expected_html = out_dir / f"{prefix}.html"

        cmd_parts = [
            "python3",
            str(worker),
            "--in",
            f'"{input_path}"',
            "--out",
            f'"{out_dir}"',
            "--prefix",
            f'"{prefix}"',
            "--comparison-mode",
            f'"{comparison_mode}"',
        ]
        cmd_str = " ".join(cmd_parts)
        print("[DEBUG] PyGenomeVizCompare worker cmd:", cmd_str)

        rc = shell.run(cmd_str, shell_mode=True)
        if rc != 0:
            files_list = self._list_dir_files(out_dir)
            raise RuntimeError(
                f"PyGenomeVizCompare worker failed with exit code {rc}.\n"
                f"Files present in {out_dir}:\n{files_list}\n"
            )

        if not expected_html.is_file():
            files_list = self._list_dir_files(out_dir)
            raise RuntimeError(
                "PyGenomeVizCompare completed but the expected HTML file was not found.\n"
                f"Expected: {expected_html}\n"
                f"Files present in {out_dir}:\n{files_list}\n"
            )

        return {"report_html": File(str(expected_html))}