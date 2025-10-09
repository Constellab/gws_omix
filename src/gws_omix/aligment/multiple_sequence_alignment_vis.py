#!/usr/bin/env python3
# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from pathlib import Path
from typing import Final

from gws_core import (
    BoolParam, ConfigParams, ConfigSpecs, File,
    InputSpec, InputSpecs, IntParam, OutputSpec, OutputSpecs,
    ResourceSet, ShellProxy, StrParam, Task, TaskInputs, TaskOutputs,
    task_decorator
)

from gws_omix.base_env.msaviz_env_task import MSAVisShellProxyHelper

@task_decorator(
    "MSAVis",
    human_name="Multiple Sequence Alignment",
    short_description="Automatic Multiple Sequence Alignment using Mafft "
)
class MSAVis(Task):
    """
    Performs multiple sequence alignment using MAFFT (if needed) and visualizes
    the resulting alignment using pyMSAviz. The visualization automatically adapts
    its layout to long alignments by splitting them into multiple pages (PNG files).

    Features:
    - Automatic MAFFT alignment (if sequences not pre-aligned)
    - Auto-detection of nucleotide vs protein color scheme
    - Pagination for long alignments
    Output:
    - <prefix>.aln.fa : aligned FASTA file
    - <prefix>.page_*.png : paginated alignment images
     """

    input_specs: Final[InputSpecs] = InputSpecs({
        "sequences_fasta": InputSpec(
            File,
            human_name="Sequences (FASTA)",
            short_description="Protein or nucleotide sequences to align/visualize"
        )
    })

    output_specs: Final[OutputSpecs] = OutputSpecs({
        "msa_outputs":   OutputSpec(ResourceSet, human_name="MSA figures (PNG)"),
        "mafft_outputs": OutputSpec(ResourceSet, human_name="MAFFT alignment & logs")
    })

    # Removed: dpi, seq_index (now fixed inside the worker)
    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(default_value="", short_description="Output prefix (default: input filename stem)"),
        "max_seqs": IntParam(default_value=500, min_value=1, short_description="Max sequences to keep"),
        "align_if_needed": BoolParam(default_value=True, short_description="Run MAFFT if input not aligned"),
        "threads": IntParam(default_value=8, min_value=1, short_description="MAFFT threads"),
    })

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        fasta: File = inputs["sequences_fasta"]

        shell: ShellProxy = MSAVisShellProxyHelper.create_proxy(self.message_dispatcher)
        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "msa_vis"
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix = (params["prefix"] or Path(fasta.path).stem).strip() or "msa"

        worker = Path(__file__).with_name("_multiple_sequence_alignment_vis.py")

        cmd = [
            "python3",
            str(worker),
            "--in",      f"\"{fasta.path}\"",
            "--out",     f"\"{out_dir}\"",
            "--prefix",  f"\"{prefix}\"",
            "--max-seqs", str(params["max_seqs"]),
            "--threads",  str(params["threads"]),
            "--align-if-needed", "1" if params["align_if_needed"] else "0",
        ]
        cmd_str = " ".join(cmd)
        print("[DEBUG] MSAVis worker cmd:", cmd_str)

        rc = shell.run(cmd_str, shell_mode=True)

        # Always build RS, even on error, so we keep logs
        msa_rs   = ResourceSet()
        mafft_rs = ResourceSet()

        # Figures (PNGs)
        for p in sorted(out_dir.glob(f"{prefix}*.png")):
            if p.is_file():
                msa_rs.add_resource(File(str(p)), p.name)

        # Alignments & logs
        for name in [f"{prefix}.aln.fa", f"{prefix}.trimmed.fa", f"{prefix}.mafft.log", f"{prefix}.log"]:
            p = out_dir / name
            if p.is_file():
                mafft_rs.add_resource(File(str(p)), p.name)

        for p in sorted(out_dir.glob(f"{prefix}.page_*.fa")):
            if p.is_file():
                mafft_rs.add_resource(File(str(p)), p.name)

        if rc:
            tail = ""
            log_path = out_dir / f"{prefix}.log"
            if log_path.is_file():
                try:
                    tail = "\n".join(log_path.read_text(errors="ignore").splitlines()[-80:])
                except Exception:
                    pass
            raise RuntimeError(f"pyMSAviz worker failed. See log tail:\n{tail}")

        # Worker said OK but no images? Bubble up context.
        if not any(out_dir.glob(f"{prefix}*.png")):
            tail = ""
            log_path = out_dir / f"{prefix}.log"
            if log_path.is_file():
                try:
                    tail = "\n".join(log_path.read_text(errors="ignore").splitlines()[-80:])
                except Exception:
                    pass
            raise RuntimeError("MSA finished but no PNGs were produced.\n" + tail)

        return {"msa_outputs": msa_rs, "mafft_outputs": mafft_rs}
