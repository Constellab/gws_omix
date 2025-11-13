#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from typing import Final
import shlex

from gws_core import (
    ConfigParams, ConfigSpecs, File,
    InputSpec, InputSpecs, OutputSpec, OutputSpecs,
    ShellProxy, StrParam, Task, TaskInputs, TaskOutputs,
    task_decorator
)
from ..dsdna_ug_to_pmol.dsdna_ug_to_pmol_env import DsdnaUgToPmolShellProxyHelper

@task_decorator(
    "bioseq_mw",
    human_name="DNA/RNA/Protein Molecular Weight Calculator",
    short_description="Compute MW for DNA/RNA (ss/ds, linear/circular, 5′ end) or protein from a sequence file."
)
class BioSeqMW(Task):
    """
    Molecular Weight (MW) Calculator — DNA, RNA, Protein
    ====================================================

    What this does
    --------------
    - DNA/RNA: Reads a nucleotide sequence (FASTA or plain), lets you pick strand (ss/ds),
    topology (linear/circular), and 5′ chemistry (hydroxyl/phosphate/triphosphate).
    Computes MW using base-specific masses + small end-corrections tuned to match vendor tools.
    For circular molecules, models ligation by subtracting H2O (condensation) at the junction.
    - Protein: Reads a protein sequence and uses ProtParam (Biopython) for average MW (with terminal H2O),
    amino-acid composition (%), aromaticity, instability index, pI, secondary structure fractions,
    and extinction coefficients.

    Core constants (average masses, Da)
    -----------------------------------
    DNA nucleotides: dA=313.209, dT=304.197, dC=289.184, dG=329.205
    RNA nucleotides:  A=329.21,  U=306.17,  C=305.18,  G=345.21
    dsDNA bp average: 618.004812 Da/bp
    dsRNA bp average: 644.574228 Da/bp
    Water (H2O):      18.01528 Da

    End-group (5′ chemistry) corrections (Da)
    -----------------------------------------
    Single-stranded DNA (per strand):
    5′-OH: −67.497 ; 5′-P: (−67.497 + 79.98) ; 5′-PPP: (−67.497 + 239.94)
    Double-stranded DNA (total per molecule):
    5′-OH: −123.92 ; 5′-P: +36.04 ; 5′-PPP: +355.96
    Single-stranded RNA (per strand):
    5′-OH: −72.10 ; 5′-P: (−72.10 + 79.98) ; 5′-PPP: (−72.10 + 159.96)

    Circular ligation rule
    ----------------------
    Closing a strand forms an extra phosphodiester bond and eliminates H2O:
    - circular ss: subtract 18.01528 Da
    - circular ds: subtract 36.03056 Da (2 × H2O)
    (Practically, circular NA behaves as 5′-phosphate at the junction.)

    Formulas (N = length in bases for ss, base pairs for ds)
    --------------------------------------------------------
    A) ssDNA:
    S = (#A×313.209 + #T×304.197 + #C×289.184 + #G×329.205)
    MW = S + Δ_end_ssDNA(5′ chemistry)
    if circular: MW -= 18.01528

    B) dsDNA (bp-average model):
    core = N_bp × 618.004812
    MW   = core + Δ_end_dsDNA(5′ chemistry)
    if circular: MW -= 36.03056

    C) ssRNA:
    S = (#A×329.21 + #U×306.17 + #C×305.18 + #G×345.21)
    MW = S + Δ_end_ssRNA(5′ chemistry)
    if circular: MW -= 18.01528

    D) dsRNA (bp-average model):
    core = N_bp × 644.574228
    MW   = core
    if circular: MW -= 36.03056
    (Note: dsRNA end-corrections are not applied in this model; bp-average captures bonding.)

    Protein MW (brief)
    ------------------
    Uses Biopython ProtParam:
    - molecular_weight(): average intact-chain mass (includes terminal H2O)
    - plus composition (% per residue)
    """


    input_specs: Final[InputSpecs] = InputSpecs({
        'input_path': InputSpec(File, human_name="Fasta File", short_description="Fasta Input File")
    })

    output_specs: Final[OutputSpecs] = OutputSpecs({
        "json": OutputSpec(File, human_name="Result JSON")
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(default_value="bioseq_mw", short_description="Output JSON prefix"),
        "type": StrParam(default_value="DNA", allowed_values=["DNA", "RNA", "PROTEIN"], short_description="Sequence type"),
        "strand": StrParam(
            default_value="",
            allowed_values=["", "ss", "ds"],
            short_description="For DNA/RNA: ss or ds. For PROTEIN: leave EMPTY ('')."
        ),
        "topology": StrParam(
            default_value="",
            allowed_values=["", "linear", "circular"],
            short_description="For DNA/RNA: linear or circular. For PROTEIN: leave EMPTY ('')."
        ),
        "five_prime": StrParam(
            default_value="",
            allowed_values=["", "hydroxyl", "phosphate", "triphosphate"],
            short_description="For DNA/RNA: 5′ end chemistry. For PROTEIN: leave EMPTY ('')."
        )
    })

    def _q(self, x) -> str: return shlex.quote(str(x))
    def _opt(self, name: str, val: str) -> str:
        v = "" if val is None else str(val)
        return f" --{name} {self._q(v)}"

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        shell: ShellProxy = DsdnaUgToPmolShellProxyHelper.create_proxy(self.message_dispatcher)

        in_file: File = inputs["input_path"]
        in_path = Path(in_file.path)
        if not in_path.is_file():
            raise FileNotFoundError(f"Input file not found: {in_path}")

        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "bioseq_mw"
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix = (params["prefix"] or "bioseq_mw").strip() or "bioseq_mw"
        worker = Path(__file__).with_name("_bioseq_mw.py")
        if not worker.exists():
            raise FileNotFoundError(f"Worker not found: {worker}")

        cmd = f"python3 {self._q(worker)} --out {self._q(out_dir)} --prefix {self._q(prefix)}"
        cmd += self._opt("type", params["type"])
        cmd += self._opt("in", str(in_path))
        cmd += self._opt("strand", params["strand"])
        cmd += self._opt("topology", params["topology"])
        cmd += self._opt("five-prime", params["five_prime"])

        print("[DEBUG] worker cmd:", cmd)
        rc = shell.run(cmd, shell_mode=True)

        json_path = out_dir / f"{prefix}.json"
        if rc and json_path.is_file():
            return {"json": File(str(json_path))}
        if rc:
            raise RuntimeError(f"bioseq_mw worker failed (rc={rc}).")
        if not json_path.is_file():
            raise FileNotFoundError(f"JSON not produced: {json_path}")

        return {"json": File(str(json_path))}
