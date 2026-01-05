#!/usr/bin/env python3
# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from pathlib import Path
from typing import Final

from gws_core import (
    BoolParam,
    ConfigParams,
    ConfigSpecs,
    File,
    InputSpec,
    InputSpecs,
    IntParam,
    OutputSpec,
    OutputSpecs,
    ResourceSet,
    ShellProxy,
    StrParam,
    Task,
    TaskInputs,
    TaskOutputs,
    task_decorator,
)

from gws_omix.base_env.msaviz_env_task import MSAVisShellProxyHelper


@task_decorator(
    "MSAVis",
    human_name="Multiple Sequence Alignment",
    short_description="Automatic Multiple Sequence Alignment using Mafft "
)
class MSAVis(Task):
    """
    Performs multiple sequence alignment using MAFFT (L-INS-i, --localpair --maxiterate 1000
    in the worker) and visualizes the resulting alignment using pyMSAviz.

    Outputs:
    - msa_outputs   : MSA figures (PNG)
    - aligned_fasta : File -> aligned FASTA principal (ex. human.aln_3.fa)
    - alignment_log : File -> log principal (ex. human_3.log)
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
        "aligned_fasta": OutputSpec(File,        human_name="Aligned FASTA (.fa)"),
        "alignment_log": OutputSpec(File,        human_name="Alignment log (.log)"),
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(
            default_value="",
            short_description="Output prefix (default: input filename stem)"
        ),
        "align_if_needed": BoolParam(
            default_value=True,
            short_description="Run MAFFT if input not aligned"
        ),
        "threads": IntParam(
            default_value=8,
            min_value=1,
            short_description="MAFFT threads (L-INS-i)"
        ),
    })

    # ------------- helpers -------------

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
        names: list[str] = sorted(p.name for p in directory.iterdir())
        if not names:
            return "<no files>"
        return "\n".join(names)

    # ------------- run -------------

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
            # max-seqs reste interne au worker (default=500), pas exposé
            "--threads",  str(params["threads"]),
            "--align-if-needed", "1" if params["align_if_needed"] else "0",
        ]
        cmd_str = " ".join(cmd)
        print("[DEBUG] MSAVis worker cmd:", cmd_str)

        rc = shell.run(cmd_str, shell_mode=True)

        msa_rs = ResourceSet()

        # 1) PNGs -> msa_outputs
        for p in sorted(out_dir.glob(f"{prefix}*.png")):
            if p.is_file():
                msa_rs.add_resource(File(str(p)), p.name)

        # 2) Détection du FASTA principal et du log principal
        aligned_fa_path = None  # type: Path | None
        align_log_path  = None  # type: Path | None

        # FASTA principal : prioriser <prefix>.aln*.fa (ex. human.aln_3.fa)
        candidates_fa = sorted(out_dir.glob(f"{prefix}.aln*.fa"))
        if not candidates_fa:
            # fallback: n'importe quel *.fa qui commence par prefix et qui n'est pas une page
            candidates_fa = [
                p for p in sorted(out_dir.glob(f"{prefix}*.fa"))
                if p.is_file() and "page_" not in p.name
            ]
        if candidates_fa:
            aligned_fa_path = candidates_fa[0]
            print(f"[INFO] Detected aligned FASTA: {aligned_fa_path}")

        # LOG principal : <prefix>*.log (ex. human_3.log)
        candidates_log = sorted(out_dir.glob(f"{prefix}*.log"))
        if candidates_log:
            align_log_path = candidates_log[0]
            print(f"[INFO] Detected alignment log: {align_log_path}")

        # log_path pour les messages d'erreur
        log_path = align_log_path or (out_dir / f"{prefix}.log")

        # 3) gestion erreurs worker / PNG manquants
        if rc:
            tail = self._read_log_tail(log_path)
            raise RuntimeError(f"pyMSAviz worker failed. See log tail:\n{tail}")

        if not any(out_dir.glob(f"{prefix}*.png")):
            tail = self._read_log_tail(log_path)
            raise RuntimeError("MSA finished but no PNGs were produced.\n" + tail)

        # 4) Vérif qu'on a bien un FASTA aligné
        if aligned_fa_path is None or not aligned_fa_path.is_file():
            files_list = self._list_dir_files(out_dir)
            tail = self._read_log_tail(log_path)
            raise RuntimeError(
                "Aligned FASTA (.fa) not found.\n"
                f"Expected something like '{prefix}.aln*.fa'.\n\n"
                f"Files present in {out_dir}:\n{files_list}\n\n"
                f"Log tail:\n{tail}"
            )

        aligned_fasta_file = File(str(aligned_fa_path))

        # 5) Log (utile mais potentiellement absent)
        if align_log_path is not None and align_log_path.is_file():
            alignment_log_file = File(str(align_log_path))
        else:
            alignment_log_file = None

        outputs: TaskOutputs = {
            "msa_outputs": msa_rs,
            "aligned_fasta": aligned_fasta_file,
            "alignment_log": alignment_log_file,
        }
        return outputs
