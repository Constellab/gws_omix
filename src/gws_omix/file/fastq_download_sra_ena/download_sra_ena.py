#!/usr/bin/env python3

# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

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

# <<< renamed helper as requested >>>
from .download_sra_ena_env import DownloadSraEnaShellProxyHelper


@task_decorator(
    "download_sra_ena",
    human_name="Download FASTQ ",
    short_description="Run `fastq-dl --accession ... --provider ena|sra --cpus N` via worker; returns FASTQs + run-info TSV + Table."
)
class FastqDLRunner(Task):


    input_specs = InputSpecs({})

    output_specs = OutputSpecs({
        "fastq_dir": OutputSpec(Folder, human_name="FASTQ folder",
                                short_description="All downloaded FASTQ(.gz)/SRA files."),
        "run_info_tbl": OutputSpec(Table, human_name="FASTQ run metadata (TSV)",
                                   short_description="Normalized ENA/SRA run info."),
    })

    config_specs: ConfigSpecs = {
        "accession": StrParam(default_value="", short_description="Accession (PRJNA..., SRP..., SRX..., SRR...)"),
        "provider":  StrParam(default_value="sra", allowed_values=["sra","ena"],short_description="Provider: ena or sra"),
        "cpus":      IntParam(default_value=4, min_value=1, short_description="CPUs for fastq-dl (esp. for SRA)"),
    }

    # Worker file (below)
    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_download_sra_ena.py"
    )

    shell_proxy: ShellProxy = None

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """Run the task."""
        accession = (params["accession"] or "").strip()
        if not accession:
            raise ValueError("Please set a valid accession (PRJNA..., SRP..., SRX..., SRR...).")

        provider = (params["provider"] or "ena").strip().lower()
        if provider not in ("ena", "sra"):
            provider = "ena"

        cpus = int(params["cpus"] or 1)
        cpus = max(cpus, 1)

        # Activate env and get a working directory
        shell_proxy = DownloadSraEnaShellProxyHelper.create_proxy(self.message_dispatcher)

        # Create output dirs inside the working dir
        base_dir = shell_proxy.working_dir
        fastq_out_dir = os.path.join(base_dir, "fastq_dl_out")
        summary_dir   = os.path.join(base_dir, "fastq_dl_summary")
        os.makedirs(fastq_out_dir, exist_ok=True)
        os.makedirs(summary_dir, exist_ok=True)

        # Call worker with 5 args (accession provider cpus out_fastq_dir summary_dir)
        cmd = f"python3 {self.python_file_path} {accession} {provider} {cpus} {fastq_out_dir} {summary_dir}"
        shell_proxy.run(cmd, shell_mode=True)
        self.shell_proxy = shell_proxy

        # Prepare outputs
        fastq_folder = Folder(fastq_out_dir)

        # Import TSV into a Table
        summary_tsv_path = os.path.join(summary_dir, "fastq-run-info.tsv")
        if not os.path.exists(summary_tsv_path):
            # Safety: create header-only TSV if worker didn't produce one (it normally does)
            with open(summary_tsv_path, "w", encoding="utf-8") as f:
                f.write("study_accession\tsample_accession\texperiment_accession\trun_accession\tscientific_name\tlibrary_layout\tfastq_md5\tfastq_ftp\tbase_count\tread_count\n")

        run_info_table = TableImporter.call(
            File(summary_tsv_path),
            {'delimiter': 'tab', 'header': 0, 'file_format': 'tsv'}
        )

        return {
            "fastq_dir": fastq_folder,
            "run_info_tbl": run_info_table,
        }

    def run_after_task(self):
        if self.shell_proxy:
            self.shell_proxy.clean_working_dir()
