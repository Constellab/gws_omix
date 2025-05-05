# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import re
import pandas as pd

from gws_core import (
    ConfigParams, Folder, Task, InputSpecs, OutputSpecs,
    TaskInputs, TaskOutputs, task_decorator,
    InputSpec, OutputSpec, ConfigSpecs, File
)
from .salmon_env import SalmonShellProxyHelper  # A shell proxy for Salmon

@task_decorator(
    "Salmon_MergeMatrix",
    human_name="Salmon_MergeMatrix",
    short_description="Combine multiple 'quant.sf' files into a gene-level raw-count matrix using a GTF annotation."
)
class SalmonMergeMatrix(Task):
    """
    This task:
      1) Takes as input the output folder from Salmon_Quant (one subfolder per sample with 'quant.sf'),
         and a GTF annotation file.
      2) Parses the GTF to build a transcript→gene mapping.
      3) For each sample, reads 'Name' & 'NumReads' from quant.sf, joins to gene mapping,
         and sums counts per gene.
      4) Merges all gene-level counts across samples (outer join).
      5) Rounds counts to integers and writes a CSV with 'gene_id' as first column.
    """
    input_specs = InputSpecs({
        'salmon_quant_folder': InputSpec(
            Folder,
            human_name="Salmon Quant Folder",
            short_description="Folder produced by Salmon_Quant, containing sample subfolders with quant.sf"
        ),
        'annotation_gtf': InputSpec(
            File,
            human_name="GTF Annotation",
            short_description="GTF file used to extract transcript-to-gene mapping"
        )
    })
    output_specs = OutputSpecs({
        'matrix': OutputSpec(
            File,
            human_name="Merged Gene-Level Counts",
            short_description="CSV matrix (gene_id vs. raw integer counts per sample)"
        )
    })
    config_specs: ConfigSpecs = {}

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        salmon_folder: Folder = inputs['salmon_quant_folder']
        gtf_file: File = inputs['annotation_gtf']

        # 1) Parse GTF to extract transcript_id → gene_id mapping
        tx2gene = {}
        pattern_tx = re.compile(r'transcript_id "([^"]+)"')
        pattern_gn = re.compile(r'gene_id "([^"]+)"')
        with open(gtf_file.path, 'r') as gf:
            for line in gf:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                attrs = fields[8]
                m_tx = pattern_tx.search(attrs)
                m_gn = pattern_gn.search(attrs)
                if m_tx and m_gn:
                    tx2gene[m_tx.group(1)] = m_gn.group(1)
        if not tx2gene:
            raise Exception("No transcript_id→gene_id mapping found in GTF.")
        tx2gene_df = pd.DataFrame(
            list(tx2gene.items()),
            columns=['transcript_id', 'gene_id']
        )

        # Prepare work directory
        shell = SalmonShellProxyHelper.create_proxy(self.message_dispatcher)
        work_dir = os.path.join(shell.working_dir, 'merged')
        os.makedirs(work_dir, exist_ok=True)

        # 2) Locate quant.sf in each sample subfolder
        samples = []
        for name in os.listdir(salmon_folder.path):
            subdir = os.path.join(salmon_folder.path, name)
            qsf = os.path.join(subdir, 'quant.sf')
            if os.path.isfile(qsf):
                samples.append((name, qsf))
        if not samples:
            raise Exception("No quant.sf files found in salmon_quant_folder.")

        # 3) Read and collapse each sample to gene-level
        gene_dfs = []
        for sample_name, qsf_path in samples:
            df = pd.read_csv(qsf_path, sep='\t', usecols=['Name', 'NumReads'], dtype={'Name': str})
            df.rename(columns={'Name': 'transcript_id', 'NumReads': sample_name}, inplace=True)
            # map transcripts to genes and sum
            merged = pd.merge(df, tx2gene_df, on='transcript_id', how='inner')
            gene_counts = merged.groupby('gene_id', as_index=False)[sample_name].sum()
            gene_dfs.append(gene_counts)

        # 4) Merge gene-level counts across all samples
        merged_counts = gene_dfs[0]
        for gdf in gene_dfs[1:]:
            merged_counts = pd.merge(merged_counts, gdf, on='gene_id', how='outer')

        # 5) Fill missing with zeros, round, convert to int
        for col in merged_counts.columns:
            if col != 'gene_id':
                merged_counts[col] = merged_counts[col].fillna(0).round().astype(int)

        # 6) Write final CSV
        out_csv = os.path.join(work_dir, 'salmon_merged_matrix.csv')
        merged_counts.to_csv(out_csv, index=False)

        return {'matrix': File(out_csv)}
