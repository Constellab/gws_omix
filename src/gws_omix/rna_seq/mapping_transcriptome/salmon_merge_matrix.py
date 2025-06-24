# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

"""
Combine multiple *quant.sf* files produced by Salmon into a gene‑level raw‑count
matrix, taking care of transcript‑version suffixes (".1", ".2", …) that may be
present in the quant.sf but absent in the GTF.

* Strip the version suffix from **both** transcript_id and gene_id when parsing
  the GTF.
* Strip the version suffix from the *Name* column of *quant.sf* before the join.
"""

import os
import re
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

    config_specs: ConfigSpecs = ConfigSpecs({})

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        salmon_folder: Folder = inputs['salmon_quant_folder']
        gtf_file: File = inputs['annotation_gtf']

        # 1) Build transcript→gene map from GTF
        tx2gene: dict[str, str] = {}
        patt_tx = re.compile(r'transcript_id "([^\"]+)"')
        patt_gn = re.compile(r'gene_id "([^\"]+)"')

        with open(gtf_file.path) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 9:
                    continue
                attrs = fields[8]
                m_tx = patt_tx.search(attrs)
                m_gn = patt_gn.search(attrs)
                if m_tx and m_gn:
                    tx_id = m_tx.group(1).split('.')[0]  # keep part before first dot
                    gn_id = m_gn.group(1).split('.')[0]
                    tx2gene[tx_id] = gn_id

        if not tx2gene:
            raise Exception("No transcript_id→gene_id pairs found in GTF. Check annotation file.")

        tx2gene_df = pd.DataFrame(tx2gene.items(), columns=['transcript_id', 'gene_id'])

        # 2) Prepare working directory inside Salmon shell environment
        shell = SalmonShellProxyHelper.create_proxy(self.message_dispatcher)
        work_dir = os.path.join(shell.working_dir, 'merged')
        os.makedirs(work_dir, exist_ok=True)

        # 3) Locate quant.sf files (one per sample)
        samples: list[tuple[str, str]] = []
        for entry in os.listdir(salmon_folder.path):
            qsf_path = os.path.join(salmon_folder.path, entry, 'quant.sf')
            if os.path.isfile(qsf_path):
                samples.append((entry, qsf_path))

        if not samples:
            raise Exception("No quant.sf files found in the provided folder.")

        # 4) Convert each sample to gene-level counts
        gene_tables: list[pd.DataFrame] = []
        for sample_name, qsf_path in samples:
            df = pd.read_csv(qsf_path, sep='\t', usecols=['Name', 'NumReads'], dtype={'Name': str})
            # Drop version suffix from transcript IDs
                        # Strip version suffix from transcript IDs (first dot only)
            df['transcript_id'] = df['Name'].str.split('.', n=1, expand=False).str[0]
            df.rename(columns={'NumReads': sample_name}, inplace=True)
            df = df[['transcript_id', sample_name]]

            merged = pd.merge(df, tx2gene_df, on='transcript_id', how='inner')
            gene_counts = merged.groupby('gene_id', as_index=False)[sample_name].sum()
            gene_tables.append(gene_counts)

        if not gene_tables:
            raise Exception("Mapping transcripts to genes produced empty tables. Check GTF vs index.")

        # 5) Merge across samples (outer join to keep all genes)
        combined = gene_tables[0]
        for tbl in gene_tables[1:]:
            combined = pd.merge(combined, tbl, on='gene_id', how='outer')


        # 6) Fill missing with 0, round, cast to int

        for col in combined.columns[1:]:  # skip gene_id column
            combined[col] = combined[col].fillna(0).round().astype(int)

        # 7) Save result
        out_path = os.path.join(work_dir, 'salmon_merged_matrix.csv')
        combined.to_csv(out_path, index=False)

        return {'matrix': File(out_path)}
