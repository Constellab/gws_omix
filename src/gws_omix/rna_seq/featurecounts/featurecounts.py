# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import re
import pandas as pd

from gws_core import (
    ConfigParams, Folder, IntParam, StrParam, Task, InputSpecs, OutputSpecs,
    TaskInputs, TaskOutputs, task_decorator, InputSpec, OutputSpec, ConfigSpecs, File
)
from .featurecounts_env import FeatureCountsShellProxyHelper



@task_decorator("FeatureCounts", human_name="FeatureCounts",
                short_description="Quantify read counts using featureCounts, then produce a final CSV with gene_id and gene_name.")
class FeatureCounts(Task):
    """
    This task runs featureCounts to count the number of reads (raw counts).
    Then it parses the raw output to create a final CSV containing:
      - A first column 'gene_id' (replacing the default 'Geneid'),
      - A second column 'gene_name' (from the GTF),
      - And columns for each sample's read counts
    """

    input_specs: InputSpecs = InputSpecs({
        'annotation_file': InputSpec(File, human_name="GTF Annotation",
                                     short_description="Reference GTF for counting"),
        'bam_files': InputSpec(Folder, human_name="BAM Files Folder",
                               short_description="Folder containing one or more BAM files")
    })

    output_specs: OutputSpecs = OutputSpecs({
        'output': OutputSpec(
            File,
            human_name="Counts Matrix",
            short_description="CSV matrix of raw counts, including gene_id and gene_name"
        )
    })

    config_specs: ConfigSpecs = ConfigSpecs({
        "threads": IntParam(default_value=4, min_value=1, short_description="Number of threads"),
        "sequencing_type": StrParam(
            allowed_values=["Paired-end", "Single-end"],
            short_description="Library type for featureCounts"
        ),
        "strandedness": IntParam(
            default_value=0,
            min_value=0,
            max_value=2,
            short_description=(
                "indicating if strand-specific read counting should be performed."
                "It has three possible values: 0 (unstranded), 1 (stranded) and 2 (reversely stranded). 0 by default."
            )
        ),
    })

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        # 1) Retrieve inputs
        annotation_file: File = inputs['annotation_file']
        bam_folder: Folder = inputs['bam_files']

        # 2) Retrieve config
        threads = params["threads"]
        seq_type = params["sequencing_type"]
        strandedness  = params["strandedness"]

        # 3) Create working directory
        shell_proxy = FeatureCountsShellProxyHelper.create_proxy(self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'featurecounts_result')
        os.makedirs(result_path, exist_ok=True)

        # 4) Collect .bam files
        bam_list = [os.path.join(bam_folder.path, f)
                    for f in os.listdir(bam_folder.path) if f.endswith('.bam')]
        if not bam_list:
            raise Exception("No BAM files found in the provided folder.")

        raw_counts_txt = os.path.join(result_path, "featurecounts_raw.txt")

        # 5) Build featureCounts command
        paired_option = "-p" if seq_type == "Paired-end" else ""
        bam_files_str = " ".join(bam_list)
        featurecounts_cmd = (
            f"featureCounts {paired_option} "
            f"-F GTF "                              # treat file as GTF
            f"-O -M --fraction "
            f"-s {strandedness} "
            f"-t exon "                             # feature type to count
            f"-g gene_id "                          # group by gene_id
            f"-T {threads} "
            f"-a {annotation_file.path} "
            f"-o {raw_counts_txt} "
            f"{bam_files_str}"
        )
        print("[DEBUG] Running command:", featurecounts_cmd)

        # 6) Run featureCounts
        res = shell_proxy.run(featurecounts_cmd, shell_mode=True)
        if res != 0:
            raise Exception("featureCounts failed")

        # 7) Build the final CSV (with gene_id, gene_name, and sample columns)
        final_csv = os.path.join(result_path, "counts_matrix_with_name.csv")
        geneid_to_name = self._parse_gtf_for_gene_names(annotation_file.path)
        self._generate_clean_matrix(raw_counts_txt, final_csv, geneid_to_name)

        return {'output': File(final_csv)}

    def _generate_clean_matrix(self, featurecounts_file: str, output_file: str,
                               geneid_to_name: dict):
        # 1) read featureCounts raw output
        df = pd.read_csv(featurecounts_file, sep='\t', comment='#')
        # 2) set 'Geneid' as index
        df.set_index('Geneid', inplace=True)

        # 3) select only sample count columns (skip meta-columns)
        count_df = df.iloc[:, 5:]

        # 4) rename columns => remove .bam and _trimmed suffix
        def clean_col(path: str) -> str:
            basename = os.path.basename(path)
            name = re.sub(r'_trimmed\.bam$', '', basename)
            return re.sub(r'\.bam$', '', name)

        new_cols = [clean_col(c) for c in count_df.columns]
        count_df.columns = new_cols

        # 5) insert 'gene_name' column
        count_df.insert(0, 'gene_name',
                        count_df.index.map(lambda gid: geneid_to_name.get(gid, gid)))

        # 6) reset index and rename
        count_df.reset_index(inplace=True)
        count_df.rename(columns={'Geneid': 'gene_id'}, inplace=True)

        # 7) reorder and write
        cols_order = ['gene_id', 'gene_name'] + new_cols
        count_df[cols_order].to_csv(output_file, index=False)
        print(f"[INFO] Created final CSV => {output_file}")

    def _parse_gtf_for_gene_names(self, gtf_path: str) -> dict:
        geneid2name = {}
        with open(gtf_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                attrs = fields[8]
                gid = self._extract_gtf_attribute(attrs, 'gene_id')
                gname = self._extract_gtf_attribute(attrs, 'gene_name')
                if gid:
                    geneid2name[gid] = gname or gid
        return geneid2name

    def _extract_gtf_attribute(self, attr_str: str, key: str) -> str:
        match = re.search(rf'{key}\s+"([^"]+)"', attr_str)
        return match.group(1) if match else None
