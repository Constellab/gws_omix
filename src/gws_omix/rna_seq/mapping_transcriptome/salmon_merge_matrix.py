# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import pandas as pd

from gws_core import (
    ConfigParams, Folder, Task, InputSpecs, OutputSpecs,
    TaskInputs, TaskOutputs, task_decorator,
    InputSpec, OutputSpec, ConfigSpecs, File
)
from .salmon_env import SalmonShellProxyHelper


@task_decorator(
    "Salmon_MergeMatrix",
    human_name="Salmon_MergeMatrix",
    short_description="Combine multiple 'quant.sf' files into a transcript-level raw-count matrix (Name -> transcript_id)."
)
class SalmonMergeMatrix(Task):
    """
    This task:
      1) Reads Salmon quant.sf files (one subfolder per sample).
      2) Builds a transcript-level raw count matrix using Salmon 'Name' as transcript IDs.
      3) Renames Salmon 'Name' into 'transcript_id' in output.
      4) Merges all samples (outer join), fills missing with 0, casts to int, writes CSV.
    """

    input_specs = InputSpecs({
        "salmon_quant_folder": InputSpec(
            Folder,
            human_name="Salmon Quant Folder",
            short_description="Folder produced by Salmon_Quant, containing sample subfolders with quant.sf"
        ),
    })

    output_specs = OutputSpecs({
        "matrix": OutputSpec(
            File,
            human_name="Merged Transcript-Level Counts",
            short_description="CSV matrix (transcript_id vs. raw integer counts per sample)"
        )
    })

    config_specs: ConfigSpecs = ConfigSpecs({})

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        salmon_folder: Folder = inputs["salmon_quant_folder"]

        # Prepare working directory inside Salmon shell environment
        shell = SalmonShellProxyHelper.create_proxy(self.message_dispatcher)
        work_dir = os.path.join(shell.working_dir, "merged")
        os.makedirs(work_dir, exist_ok=True)

        # Locate quant.sf files
        samples: list[tuple[str, str]] = []
        for entry in os.listdir(salmon_folder.path):
            qsf_path = os.path.join(salmon_folder.path, entry, "quant.sf")
            if os.path.isfile(qsf_path):
                samples.append((entry, qsf_path))

        if not samples:
            raise Exception("No quant.sf files found in the provided folder.")

        # Build per-sample transcript tables
        tx_tables: list[pd.DataFrame] = []
        for sample_name, qsf_path in samples:
            df = pd.read_csv(
                qsf_path,
                sep="\t",
                usecols=["Name", "NumReads"],
                dtype={"Name": str},
            )
            df.rename(columns={"Name": "transcript_id", "NumReads": sample_name}, inplace=True)

            # Sum duplicates if any
            df = df.groupby("transcript_id", as_index=False)[sample_name].sum()
            tx_tables.append(df)

        # Merge across samples (outer join)
        combined = tx_tables[0]
        for tbl in tx_tables[1:]:
            combined = pd.merge(combined, tbl, on="transcript_id", how="outer")

        # Fill missing with 0 and cast to int
        for col in combined.columns[1:]:
            combined[col] = combined[col].fillna(0).round().astype(int)

        # Save
        out_path = os.path.join(work_dir, "salmon_merged_matrix.csv")
        combined.to_csv(out_path, index=False)

        return {"matrix": File(out_path)}
