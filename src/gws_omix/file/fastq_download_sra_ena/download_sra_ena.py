#!/usr/bin/env python3

# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import shlex

from gws_core import (
    ConfigParams,
    ConfigSpecs,
    File,
    Folder,
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

from .download_sra_ena_env import DownloadSraEnaShellProxyHelper


@task_decorator(
    "download_sra_ena",
    human_name="Download FASTQ ",
    short_description="Download FASTQ files from SRA or ENA repositories."
)
class FastqDLRunner(Task):
    """
    Python wrapper around `fastq-dl` (https://github.com/rpetit3/fastq-dl) that downloads FASTQ files for one or more accessions
    (ENA or SRA), then merges/creates a `fastq-run-info.tsv` summary and keeps only key columns.

    Inputs:
        accession (str): One or more accession IDs separated by commas
                         (e.g. "SRX1,SRX2,PRJNA...,SRR...")
        provider (str): "ena" or "sra" (default/fallback: "ena")
        cpus (int): number of CPU threads (default: 2)

    Outputs:
        - Directory containing downloaded FASTQ file(s) (e.g. *.fastq.gz)
          (all saved in the same folder)
        - `fastq-run-info.tsv` table listing all runs resolved under all input accessions
    """

    input_specs = InputSpecs({})

    output_specs = OutputSpecs({
        "fastq_dir": OutputSpec(
            Folder,
            human_name="FASTQ folder",
            short_description="All downloaded FASTQ(.gz)/SRA files."
        ),
        "run_info_tbl": OutputSpec(
            Table,
            human_name="FASTQ run metadata (TSV)",
            short_description="Normalized ENA/SRA run info."
        ),
    })

    config_specs: ConfigSpecs = {
        "accession": StrParam(
            default_value="",
            short_description="Accession(s) comma-separated (PRJNA..., SRP..., SRX..., SRR...)"
        ),
        "provider": StrParam(
            default_value="sra",
            allowed_values=["sra", "ena"],
            short_description="Provider: ena or sra"
        ),
        "cpus": IntParam(
            default_value=2,
            min_value=1,
            short_description="CPUs for fastq-dl (esp. for SRA)"
        ),
    }

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_download_sra_ena.py"
    )

    shell_proxy: ShellProxy = None

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        accession_raw = (params["accession"] or "").strip()
        if not accession_raw:
            raise ValueError("Please set at least one valid accession (comma-separated).")

        # Normalize comma-separated accessions (remove blanks/spaces)
        accessions = [a.strip() for a in accession_raw.split(",") if a.strip()]
        if not accessions:
            raise ValueError("Please set at least one valid accession (comma-separated).")

        accession_arg = ",".join(accessions)

        provider = (params["provider"] or "ena").strip().lower()
        if provider not in ("ena", "sra"):
            provider = "ena"

        cpus = int(params["cpus"] or 1)
        cpus = max(cpus, 1)

        shell_proxy = DownloadSraEnaShellProxyHelper.create_proxy(self.message_dispatcher)

        # Create output dirs inside the working dir
        base_dir = shell_proxy.working_dir
        fastq_out_dir = os.path.join(base_dir, "fastq_dl_out")
        summary_dir = os.path.join(base_dir, "fastq_dl_summary")
        os.makedirs(fastq_out_dir, exist_ok=True)
        os.makedirs(summary_dir, exist_ok=True)

        # Build command safely (handles spaces in paths / args)
        cmd = "python3 {py} {acc} {prov} {cpus} {outd} {sumd}".format(
            py=shlex.quote(self.python_file_path),
            acc=shlex.quote(accession_arg),
            prov=shlex.quote(provider),
            cpus=shlex.quote(str(cpus)),
            outd=shlex.quote(fastq_out_dir),
            sumd=shlex.quote(summary_dir),
        )

        shell_proxy.run(cmd, shell_mode=True)
        self.shell_proxy = shell_proxy

        # Prepare outputs
        fastq_folder = Folder(fastq_out_dir)

        # Import merged TSV into a Table
        summary_tsv_path = os.path.join(summary_dir, "fastq-run-info.tsv")
        if not os.path.exists(summary_tsv_path):
            # Safety: create header-only TSV if worker didn't produce one (it normally does)
            with open(summary_tsv_path, "w", encoding="utf-8") as f:
                f.write(
                    "study_accession\tsample_accession\texperiment_accession\trun_accession\t"
                    "scientific_name\tlibrary_layout\tfastq_md5\tfastq_ftp\tbase_count\tread_count\n"
                )

        run_info_table = TableImporter.call(
            File(summary_tsv_path),
            {"delimiter": "tab", "header": 0, "file_format": "tsv"},
        )

        return {
            "fastq_dir": fastq_folder,
            "run_info_tbl": run_info_table,
        }

    def run_after_task(self):
        if self.shell_proxy:
            self.shell_proxy.clean_working_dir()
