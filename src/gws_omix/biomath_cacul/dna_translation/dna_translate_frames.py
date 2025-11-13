#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from typing import Final

from gws_core import (
    ConfigParams, ConfigSpecs, File,
    InputSpec, InputSpecs, OutputSpec, OutputSpecs,
    ShellProxy, StrParam, Task, TaskInputs, TaskOutputs,
    ResourceSet, task_decorator, TableImporter, Table
)

from ..dsdna_ug_to_pmol.dsdna_ug_to_pmol_env import DsdnaUgToPmolShellProxyHelper


@task_decorator(
    "dna_translate_frames_orf",
    human_name="DNA translation",
    short_description="Translate forward frames (1–3) and list ORFs (ATG→stop) per frame; outputs a text report + per-frame tables."
)
class DNATranslateFramesORF(Task):
    """
    DNA → Protein (Forward Frames 1–3) + ORFs
    =========================================
    What:
    • Translates the input DNA sequence in forward frames 1, 2, and 3 using the Standard genetic code.
    • Detects ORFs in each frame as: ATG (start) → nearest in-frame stop (TAA/TAG/TGA).
    • Writes a human-readable report with only the raw frame translations (no ORF lists),
        and three per-frame TSV tables with ORFs (header on the first line).

    Why:
    • Quickly scan all forward reading frames.
    • See high-level protein translations and structured ORF calls for downstream filtering.

    Input:
    • A DNA sequence file (FASTA or plain text). Only the first sequence is read.
    • Characters: IUPAC DNA letters (A, C, G, T, plus N and ambiguity codes). Non-letters/whitespace ignored.

    Output:
    • <prefix>.txt
        – Raw protein translations for frames 1, 2, 3 (stop codons as '*').
    • <prefix>.frame{1,2,3}_orfs.tsv (per frame)
        – Columns: frame, start_nt, end_nt, len_nt, len_aa, stop, peptide

    ORF definition:
    • A run that starts at an ATG (start codon) and ends at the first in-frame stop (TAA/TAG/TGA), inclusive.
    • Peptide is translated from start to just before the stop codon.

    Notes:
    • Only forward strand is considered.
    • Standard genetic code is used (Biopython translate, no CDS enforcement).
    """



    input_specs: Final[InputSpecs] = InputSpecs({
        "input_path": InputSpec(
            File,
            human_name="DNA Sequence File",
            short_description="FASTA or plain DNA sequence. Only the first sequence is read."
        )
    })

    output_specs: Final[OutputSpecs] = OutputSpecs({
        "report": OutputSpec(File, human_name="Text Report"),
        "orf_tables": OutputSpec(ResourceSet, human_name="Per-frame ORF Tables")
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(default_value="dna_translate_orf", short_description="Output prefix (no JSON)")
    })

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        shell: ShellProxy = DsdnaUgToPmolShellProxyHelper.create_proxy(self.message_dispatcher)

        in_file: File = inputs["input_path"]
        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "dna_translate_frames_orf"
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix = (params["prefix"] or "dna_translate_orf").strip() or "dna_translate_orf"

        worker = Path(__file__).with_name("_dna_translate_frames.py")
        if not worker.exists():
            raise FileNotFoundError(f"Worker not found: {worker}")

        cmd = f'python3 "{worker}" --in "{in_file.path}" --out "{out_dir}" --prefix "{prefix}"'
        print("[DEBUG] worker cmd:", cmd)
        rc = shell.run(cmd, shell_mode=True)

        report_path = out_dir / f"{prefix}.txt"
        tsv_paths = [
            out_dir / f"{prefix}.frame1_orfs.tsv",
            out_dir / f"{prefix}.frame2_orfs.tsv",
            out_dir / f"{prefix}.frame3_orfs.tsv",
        ]

        rs = ResourceSet()
        for p in tsv_paths:
            if p.is_file():
                try:
                    tbl: Table = TableImporter.call(
                        File(str(p)),
                        {
                            "delimiter": "tab",
                            "header": 0,        # <-- Header is on FIRST line in our TSVs
                            "file_format": "tsv",
                            "index_column": None,
                            "comment": "#"
                        }
                    )
                    rs.add_resource(tbl, p.name)
                except Exception as e:
                    print(f"[WARN] Table import failed for {p.name}: {e}. Falling back to File.")
                    rs.add_resource(File(str(p)), p.name)

        if rc:
            if report_path.is_file():
                return {"report": File(str(report_path)), "orf_tables": rs}
            raise RuntimeError(f"dna_translate_frames_orf worker failed (rc={rc}).")

        if not report_path.is_file():
            raise FileNotFoundError(f"Report not produced: {report_path}")

        return {"report": File(str(report_path)), "orf_tables": rs}
