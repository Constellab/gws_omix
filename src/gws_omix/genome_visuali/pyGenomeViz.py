#!/usr/bin/env python3
# LICENSE
# Exclusive property of Gencovery SAS
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# Written by Gencovery SAS, 2024

from pathlib import Path
from typing import Final, List, Union

from gws_core import (
    ConfigParams, ConfigSpecs, File, Folder,
    InputSpec, InputSpecs, OutputSpec, OutputSpecs,
    ShellProxy, StrParam, IntParam, BoolParam, Task, TaskInputs, TaskOutputs,
    task_decorator
)

from .pygenomeviz_env import PyGenomeVizShellProxyHelper


@task_decorator(
    "PyGenomeViz",
    human_name="Linear Genome Visualization (pyGenomeViz)",
    short_description="Interactive visualization of multiple genome files (GenBank or GFF)",hide=True
)
class PyGenomeViz(Task):
    """
    Visualizes one or more genome files using pyGenomeViz core API.

    **Input**
    - `genome_input` (File or Folder): A genome file (GenBank or GFF) or a folder containing genome files
      - GenBank: .gb, .gbk, .genbank
      - GFF: .gff, .gff3

    **Output**
    - `interactive_html` (File): Interactive HTML visualization with zoom/pan controls

    **Interactive Features**
    - **Zoom**: Mouse wheel to zoom in/out on regions
    - **Pan**: Click and drag to navigate along the genome
    - **Reset**: Double-click to reset to full view
    - **Export**: Download as PNG from the toolbar
    - **Color picker**: Click pen button to select color, then click on features to colorize them

    **Genome Comparison**
    When multiple genomes are provided, enable comparison to visualize sequence similarities:
    - **MUMmer**: Fast nucleotide alignment (recommended for closely related genomes)
    - **BLAST**: Comprehensive similarity search (works for any sequence type)
    - Links between similar regions will be displayed automatically
    """

    # ---------- I/O specs ----------
    input_specs: Final[InputSpecs] = InputSpecs({
        "genome_input": InputSpec(
            (Folder, File),
            human_name="Genome input",
            short_description="A genome file (GenBank/GFF) or a folder containing one or more genome files"
        )
    })

    output_specs: Final[OutputSpecs] = OutputSpecs({
        "interactive_html": OutputSpec(
            File,
            human_name="Interactive HTML output"
        )
    })

    # ---------- configuration ----------
    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(
            default_value="pygenomeviz",
            short_description="Output HTML file prefix"
        ),
        "enable_comparison": BoolParam(
            default_value=False,
            short_description="Enable genome comparison (requires 2+ genomes)",
            human_name="Enable comparison"
        ),
        "comparison_method": StrParam(
            default_value="mummer",
            allowed_values=["mummer", "blast"],
            short_description="Comparison method: mummer (fast) or blast (comprehensive)",
            human_name="Comparison method"
        ),
        "min_alignment_length": IntParam(
            default_value=100,
            min_value=1,
            short_description="Minimum alignment length to display comparison links",
            human_name="Min alignment length"
        )
    })

    # ---------- helpers ----------
    def _read_log_tail(self, log_path: Path, max_lines: int = 80) -> str:
        if not log_path.is_file():
            return ""
        try:
            return "\n".join(log_path.read_text(errors="ignore").splitlines()[-max_lines:])
        except Exception:
            return ""

    def _list_dir_files(self, directory: Path) -> str:
        if not directory.is_dir():
            return "<directory does not exist>"
        names: List[str] = sorted(p.name for p in directory.iterdir())
        if not names:
            return "<no files>"
        return "\n".join(names)

    # ---------- main run ----------
    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        genome_input: Union[Folder, File] = inputs["genome_input"]

        # Create shell proxy
        shell: ShellProxy = PyGenomeVizShellProxyHelper.create_proxy(self.message_dispatcher)

        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "pygenomeviz"
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix = (params["prefix"] or "pygenomeviz").strip() or "pygenomeviz"
        worker = Path(__file__).with_name("_pyGenomeViz.py")
        expected_html = out_dir / f"{prefix}.html"

        input_path = genome_input.path

        # Determine comparison mode
        compare_mode = "none"
        if params.get("enable_comparison", False):
            compare_mode = params.get("comparison_method", "mummer")

        min_length = params.get("min_alignment_length", 100)

        cmd_parts = [
            "python3",
            str(worker),
            "--in", f"\"{input_path}\"",
            "--out", f"\"{out_dir}\"",
            "--prefix", f"\"{prefix}\"",
            "--compare", compare_mode,
            "--min-length", str(min_length),
        ]
        cmd_str = " ".join(cmd_parts)
        print("[DEBUG] pyGenomeViz worker cmd:", cmd_str)

        # Run worker
        rc = shell.run(cmd_str, shell_mode=True)

        # Define log path
        wrapper_log_path = out_dir / f"{prefix}.pygenomeviz_wrapper.log"

        # Error handling
        if rc != 0:
            tail = self._read_log_tail(wrapper_log_path)
            files_list = self._list_dir_files(out_dir)
            raise RuntimeError(
                f"pyGenomeViz worker failed with exit code {rc}.\n"
                f"Files present in {out_dir}:\n{files_list}\n"
                f"Log tail:\n{tail}"
            )

        # Ensure HTML file exists
        if not expected_html.is_file():
            tail = self._read_log_tail(wrapper_log_path)
            files_list = self._list_dir_files(out_dir)
            raise RuntimeError(
                f"pyGenomeViz completed without generating expected HTML file: {expected_html.name}\n"
                f"Files present in {out_dir}:\n{files_list}\n"
                f"Log tail:\n{tail}"
            )

        return {
            "interactive_html": File(str(expected_html))
        }