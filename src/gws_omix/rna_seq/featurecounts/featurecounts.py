# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import re
import pandas as pd  # <-- We use pandas to parse and transform the featureCounts output

from gws_core import (
    ConfigParams, Folder, IntParam, StrParam, Task, InputSpecs, OutputSpecs,
    TaskInputs, TaskOutputs, task_decorator, InputSpec, OutputSpec, ConfigSpecs, File
)

from .featurecounts_env import FeatureCountsShellProxyHelper

@task_decorator("FeatureCounts", human_name="FeatureCounts",
                short_description="Quantify read counts using featureCounts")
class FeatureCounts(Task):
    """
    This task runs featureCounts to count the number of reads (raw counts).
    These counts are generated from one or more BAM files mapping to genomic features (defined in a GTF file).
    They are then converted into a clean CSV count matrix for differential analysis.

    In this version, we specifically:
      - use '-F GTF' to tell featureCounts it's reading GTF
      - use '-t CDS' so it only counts lines with 'CDS' in the 3rd column
      - use '-g gene_id' to group counts by gene_id (can replace with 'transcript_id' if needed).
    """

    input_specs: InputSpecs = InputSpecs({
        'annotation_file': InputSpec(File, human_name="GTF Annotation",
                                     short_description="Reference GTF for counting"),
        'bam_files': InputSpec(Folder, human_name="BAM Files Folder",
                               short_description="Folder containing one or more BAM files")
    })

    output_specs: OutputSpecs = OutputSpecs({
        'output': OutputSpec(File, human_name="Counts Matrix",
                             short_description="Cleaned CSV matrix of raw counts")
    })

    config_specs: ConfigSpecs = {
        # For multi-threading
        "threads": IntParam(default_value=4, min_value=1, short_description="Number of threads"),
        # Single-end or Paired-end?
        "sequencing_type": StrParam(
            allowed_values=["Paired-end", "Single-end"],
            short_description="Library type for featureCounts"
        )
    }

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """Run featureCounts, then parse and convert the output into a clean CSV count matrix."""
        annotation_file: File = inputs['annotation_file']
        bam_folder: Folder = inputs['bam_files']

        threads = params["threads"]
        seq_type = params["sequencing_type"]

        # Create working directory
        shell_proxy = FeatureCountsShellProxyHelper.create_proxy(self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'featurecounts_result')
        os.makedirs(result_path, exist_ok=True)

        # Collect all .bam files from the folder
        bam_list = []
        for fname in os.listdir(bam_folder.path):
            if fname.endswith(".bam"):
                bam_list.append(os.path.join(bam_folder.path, fname))
        if not bam_list:
            raise Exception("No BAM files found in the provided folder.")

        # The raw featureCounts output file (tab-delimited)
        raw_counts_txt = os.path.join(result_path, "featurecounts_raw.txt")

        # Build featureCounts command
        paired_option = "-p" if seq_type == "Paired-end" else ""
        bam_files_str = " ".join(bam_list)
        featurecounts_cmd = (
            f"featureCounts {paired_option} "
            f"-t CDS "           # <-- count lines where column 3 == CDS
            f"-T {threads} "
            f"-a {annotation_file.path} "
            f"-o {raw_counts_txt} "
            f"{bam_files_str}"
        )

        print("[DEBUG] Running command:", featurecounts_cmd)
        res = shell_proxy.run(featurecounts_cmd, shell_mode=True)
        if res != 0:
            raise Exception("featureCounts failed")

        # Now we parse featurecounts_raw.txt and produce "counts_matrix.csv"
        counts_matrix_csv = os.path.join(result_path, "counts_matrix.csv")
        self._generate_clean_matrix(raw_counts_txt, counts_matrix_csv)

        # Return the cleaned CSV as final output
        return {'output': File(counts_matrix_csv)}

    def _generate_clean_matrix(self, featurecounts_file: str, output_file: str):
        """
        Read the multi-sample featureCounts output file (featurecounts_file),
        rename columns by removing path and "_trimmed.bam",
        then write a final CSV count matrix to output_file.
        """

        # 1) Read the featureCounts table, skipping lines starting with '#'
        df = pd.read_csv(featurecounts_file, sep='\t', comment='#')

        # 2) Use 'Geneid' as row index
        df.set_index('Geneid', inplace=True)

        # 3) The first 5 columns after 'Geneid' are (Chr, Start, End, Strand, Length).
        #    The actual raw counts start at column index 5.
        #    We'll keep only these read count columns.
        count_df = df.iloc[:, 5:]  # columns from 6th onward

        # 4) Function to transform
        #    /path/to/SRR13978641_trimmed.bam => SRR13978641
        def clean_col_name(full_bam_path: str) -> str:
            # Extract filename => SRR13978641_trimmed.bam
            basename = os.path.basename(full_bam_path)
            # Remove "_trimmed.bam" => SRR13978641
            no_trimmed = re.sub(r'_trimmed\.bam$', '', basename)
            return no_trimmed

        # 5) Clean each column header
        new_cols = []
        for col in count_df.columns:
            new_cols.append(clean_col_name(col))

        count_df.columns = new_cols

        # 6) Write out to CSV
        #    If you want a TSV, use sep='\t'.
        count_df.to_csv(output_file, sep=',')
        print(f"[INFO] Wrote cleaned count matrix to {output_file}. Columns: {new_cols}")
