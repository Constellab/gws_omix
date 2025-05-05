from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pandas as pd
import scanpy as sc
import numpy as np
import seaborn as sns
import decoupler as dc
import matplotlib.pyplot as plt
import sys

class DeseqAnalyzer:
    def __init__(
        self,
        count_table_path,
        metadata_path,
        genes_colname,
        control_condition=None,
        unnormal_condition=None,
        pvalue_value=None,
        log2FoldChange_value=None
    ):
        self.count_table_path = count_table_path
        self.metadata_path = metadata_path
        self.genes_colname = genes_colname  # can be 'gene_id' or 'gene_name'
        self.control_condition = control_condition
        self.unnormal_condition = unnormal_condition
        self.pvalue_value = pvalue_value
        self.log2FoldChange_value = log2FoldChange_value
        self.sigs = None
        self.dds = None

    def load_data(self):
        """
        1) Reads the count matrix from CSV.
        2) Sets either 'gene_id' or 'gene_name' as index (depending on self.genes_colname).
        3) Drops the other text column if it exists, so summation won't fail.
        4) Converts remaining columns to numeric, removing rows with total counts=0.
        5) Transposes so samples become rows, and aligns them with metadata.
        6) Creates DeseqDataSet.
        """
        # 1) Read the CSV (count matrix)
        counts = pd.read_csv(self.count_table_path)

        # 2) Set the user-chosen column as the DataFrame index
        #    (for example, if genes_colname='gene_name', we use that as index)
        counts.set_index(self.genes_colname, inplace=True)

        # 3) Conditionally drop the other column, so it won't be included in numeric summations
        possible_text_cols = ["gene_id", "gene_name"]
        for col in possible_text_cols:
            # Drop if present and it's NOT the user-chosen index
            if col in counts.columns and col != self.genes_colname:
                counts.drop(columns=[col], inplace=True)

        # 4) Convert all remaining columns to numeric, and drop rows whose total = 0
        counts = counts.apply(pd.to_numeric, errors='coerce').fillna(0)
        # Cast the count values to integers (round them, if necessary, to avoid decimals)
        counts = counts.astype(int)
        counts = counts[counts.sum(axis=1) > 0]

        # 5) Transpose: now rows are samples, columns are genes
        counts = counts.T

        # Read metadata, set sample column as index, align
        metadata = pd.read_csv(self.metadata_path, sep='\t', comment='#')
        metadata = metadata.set_index('Sample')
        counts = counts.reindex(metadata.index)

        # 6) Create PyDESeq2 dataset
        self.dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors="Condition")

    def run_deseq_stats(self):
        """Runs DESeq2, then applies user-defined pvalue/log2foldChange thresholds."""
        self.dds.deseq2()  # Fit the model

        # Perform Wald test using user-specified conditions
        stat_res = DeseqStats(
            self.dds,
            contrast=('Condition', self.unnormal_condition, self.control_condition)
        )
        stat_res.summary()
        res = stat_res.results_df

        # Filter based on user-supplied thresholds
        if self.pvalue_value is not None and self.log2FoldChange_value is not None:
            self.sigs = res[
                (res.pvalue < self.pvalue_value)
                & (abs(res.log2FoldChange) > self.log2FoldChange_value)
            ]
        else:
            self.sigs = res

        # Sort by log2FoldChange descending
        self.sigs = self.sigs.sort_values(by='log2FoldChange', ascending=False)

    def save_desq2_results_table(self):
        """Save the filtered DE results to CSV."""
        self.sigs.to_csv('pydesq2_results_table.csv', index=True)

    def perform_pca(self):
        """Normalize for PCA, compute PCA."""
        self.dds.layers['counts'] = self.dds.X.copy()
        sc.pp.normalize_total(self.dds)
        sc.tl.pca(self.dds)

    def save_pca_metadata(self):
        """Export the first 2 PC coordinates along with the metadata."""
        pca_coords = pd.DataFrame(self.dds.obsm['X_pca'][:, :2], columns=['PC1', 'PC2'])
        metadata = pd.read_csv(self.metadata_path, sep='\t', comment='#')

        # Merge PCA coords with metadata
        pca_metadata = pd.concat([pca_coords, metadata], axis=1, join='inner')
        pca_metadata.to_csv('pca_metadata.csv', index=False)

    def save_pca_proportions(self):
        """Save the proportion of variance explained by PC1 and PC2."""
        pca_proportions = self.dds.uns['pca']['variance_ratio'][:2]
        pca_proportions_df = pd.DataFrame({
            'PC1 Proportion': [pca_proportions[0]],
            'PC2 Proportion': [pca_proportions[1]]
        })
        pca_proportions_df.to_csv('pca_proportions.csv', index=False)

    def generate_clustermap(self):
        """Produce a Seaborn clustermap (Heatmap.png)."""
        # Normalized counts (log1p)
        self.dds.layers['log1p'] = np.log1p(self.dds.layers['normed_counts'])

        # Slice out only the 'significant' genes
        dds_sigs = self.dds[:, self.sigs.index]
        grapher = pd.DataFrame(
            dds_sigs.layers['log1p'].T,
            index=dds_sigs.var_names,
            columns=dds_sigs.obs_names
        )
        # For Plotly heatmap creation outside this script
        grapher.to_csv('grapher.csv', index=True)

        # Create a seaborn clustermap for a static PNG
        clustermap = sns.clustermap(grapher, z_score=0, cmap='RdYlBu_r')
        clustermap.savefig('Heatmap.png')

    def generate_volcano_plot(self):
        """Create a volcano plot PNG using decoupler (top=30)."""
        dc.plot_volcano_df(self.sigs, x='log2FoldChange', y='pvalue', top=30)
        plt.savefig('volcano_plot.png')

    def run_analysis(self):
        """Master method that runs all steps end-to-end."""
        self.load_data()
        self.run_deseq_stats()
        self.perform_pca()
        self.save_pca_metadata()
        self.save_pca_proportions()
        self.save_desq2_results_table()
        self.generate_clustermap()
        self.generate_volcano_plot()

# =============================
# Command-line entry point
# =============================
if __name__ == "__main__":
    count_table_file = sys.argv[1]
    metadata_file = sys.argv[2]
    genes_colname = str(sys.argv[3])  # either "gene_id" or "gene_name"
    control_condition = str(sys.argv[4])
    unnormal_condition = str(sys.argv[5])
    pvalue_value = float(sys.argv[6])
    log2FoldChange_value = float(sys.argv[7])

    deseq_analyzer = DeseqAnalyzer(
        count_table_file,
        metadata_file,
        genes_colname=genes_colname,
        control_condition=control_condition,
        unnormal_condition=unnormal_condition,
        pvalue_value=pvalue_value,
        log2FoldChange_value=log2FoldChange_value
    )

    deseq_analyzer.run_analysis()
