#!/usr/bin/env python3
import os
from pathlib import Path
from typing import Final

from gws_core import (
    ConfigParams,
    ConfigSpecs,
    File,
    InputSpec,
    InputSpecs,
    IntParam,
    OutputSpec,
    OutputSpecs,
    ShellProxy,
    StrParam,
    Table,
    TableImporter,
    Task,
    TaskInputs,
    TaskOutputs,
    task_decorator,
)

from gws_omix.base_env.eggnog_env_task import EggnogShellProxyHelper


@task_decorator(
    "Emapper",
    human_name="eggNOG Mapper",
    short_description="Functional annotation with eggNOG-mapper using DIAMOND."
)
class EggnogMapperTask(Task):
    """
   The EggNOG script is a command-line tool that automates the use of eggNOG-mapper, a high-throughput functional annotation tool
   based on orthology assignments from the eggNOG database. eggNOG-mapper allows researchers to functionally annotate protein, transcript,
   genome, or metagenome sequences with high accuracy, avoiding false positives from close paralogs by focusing on fine-grained orthologs.
   This script simplifies the process by handling all preparatory steps:
   it checks and downloads required core databases (like eggnog.db) and the DIAMOND protein database,
   ensures the environment is properly set up, and executes the emapper.py pipeline using user-defined parameters.
   The script takes four inputs—an input FASTA file, an output directory, the number of CPUs,
   and the sequence type (proteins, CDS, genome, or metagenome)—and outputs a cleaned annotation.tsv table with meaningful gene annotations.
    This automation makes eggNOG-mapper easier to integrate into large-scale genome or transcriptome annotation workflows.
     """

    input_specs: Final[InputSpecs] = InputSpecs({
        "fasta": InputSpec(
            File, human_name="FASTA file",
            short_description="Protein or transcript FASTA to annotate"),
    })

    output_specs: Final[OutputSpecs] = OutputSpecs({
        "annotation_table": OutputSpec(
            Table,
            human_name="eggNOG annotations",
            short_description="Cleaned .tsv table output from eggNOG-mapper"),
    })


    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "itype": StrParam(
            default_value="proteins",
            allowed_values=["proteins", "CDS", "genome", "metagenome"],
            short_description="Sequence type for eggNOG input (e.g., proteins, CDS, etc.)"),
        "cpus": IntParam(
            default_value=25, min_value=1,
            short_description="Number of CPU threads to use"),
    })

    python_file_path: Final[str] = os.path.join(
        os.path.abspath(os.path.dirname(__file__)), "_eggnog_emapper.py"
    )

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        fasta: File = inputs["fasta"]
        itype: str = params["itype"]
        cpus: int = params["cpus"]

        # Create shell proxy with activated eggnog environment
        shell: ShellProxy = EggnogShellProxyHelper.create_proxy(self.message_dispatcher)
        work_dir = Path(shell.working_dir)
        prefix = Path(fasta.path).stem
        output_dir = work_dir / f"{prefix}_output"

        cmd = (
            f"python3 {self.python_file_path} "
            f"{fasta.path} {output_dir} {cpus} {itype}"
        )

        ret_code = shell.run(cmd, shell_mode=True)
        if ret_code != 0:
            raise RuntimeError("eggNOG-mapper pipeline failed. Check logs for more info.")

        # Check expected output
        cleaned_output_path = output_dir / "annotation.tsv"
        if not cleaned_output_path.exists():
            raise FileNotFoundError(f"Expected cleaned annotation file not found: {cleaned_output_path}")

        annotation_table = TableImporter.call(
            File(str(cleaned_output_path)),
            {
                "delimiter": "tab",
                "header": -1,
                "file_format": "tsv",
                "index_column": None,
                "comment": "#"
            }
        )

        # Rename columns from default "0", "1", ..., "20" to meaningful headers
        column_names = [
            "Query", "seed_Ortholog", "E-value", "Score", "eggNOG_OGs", "max_annot_lvl", "COG_category",
            "Description", "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway",
            "KEGG_Module", "KEGG_Reaction", "KEGG_Rclass", "BRITE", "KEGG_TC", "CAZy",
            "BiGG_Reaction", "PFAMs"
        ]

        if len(annotation_table.column_names) != len(column_names):
            raise ValueError("Mismatch between detected and expected column count")

        for i, new_col in enumerate(column_names):
            annotation_table.set_column_name(str(i), new_col)

        return {
            "annotation_table": annotation_table
        }
