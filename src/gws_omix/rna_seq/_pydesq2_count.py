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
        self.genes_colname = genes_colname
        self.control_condition = control_condition
        self.unnormal_condition = unnormal_condition
        self.pvalue_value = pvalue_value
        self.log2FoldChange_value = log2FoldChange_value
        self.sigs = None
        self.dds = None

    def load_data(self):
        # Load count matrix and metadata
        counts = pd.read_csv(self.count_table_path)
        metadata = pd.read_csv(self.metadata_path, sep='\t', comment='#')

        # Set geneid as index
        counts = counts.set_index(self.genes_colname)

        # Remove rows with zero counts
        counts = counts[counts.sum(axis=1) > 0]

        # Transpose the matrix (samples become rows)
        counts = counts.T

        # Set the sample column as index in metadata
        metadata = metadata.set_index('Sample')

        # Reindex counts to match metadata order (optional but recommended)
        counts = counts.reindex(metadata.index)

        # Create DeseqDataSet
        self.dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors="Condition")

    def run_deseq_stats(self):
        # Run DESeq2
        self.dds.deseq2()

        # Perform Wald test
        stat_res = DeseqStats(self.dds, contrast=('Condition', self.unnormal_condition, self.control_condition))
        stat_res.summary()
        res = stat_res.results_df

        # Filter based on user-supplied pvalue and log2FoldChange thresholds
        if self.pvalue_value is not None and self.log2FoldChange_value is not None:
            self.sigs = res[
                (res.pvalue < self.pvalue_value)
                & (abs(res.log2FoldChange) > self.log2FoldChange_value)
            ]
        else:
            self.sigs = res

        # Sort by log2FoldChange in descending order (from positive to negative)
        self.sigs = self.sigs.sort_values(by='log2FoldChange', ascending=False)

    def save_desq2_results_table(self):
        # Save filtered genes to CSV
        self.sigs.to_csv('pydesq2_results_table.csv', index=True)

    def perform_pca(self):
        # Normalization for PCA
        self.dds.layers['counts'] = self.dds.X.copy()
        sc.pp.normalize_total(self.dds)
        sc.tl.pca(self.dds)

    def save_pca_metadata(self):
        # Access the PCA coordinates
        pca_coords = pd.DataFrame(self.dds.obsm['X_pca'][:, :2], columns=['PC1', 'PC2'])
        metadata = pd.read_csv(self.metadata_path, sep='\t', comment='#')
        pca_metadata = pd.concat([pca_coords, metadata], axis=1, join='inner')

        # Save PCA metadata
        pca_metadata.to_csv('pca_metadata.csv', index=False)

    def save_pca_proportions(self):
        # Access PCA proportions
        pca_proportions = self.dds.uns['pca']['variance_ratio'][:2]
        pca_proportions_df = pd.DataFrame({
            'PC1 Proportion': [pca_proportions[0]],
            'PC2 Proportion': [pca_proportions[1]]
        })
        pca_proportions_df.to_csv('pca_proportions.csv', index=False)

    def generate_clustermap(self):
        # Normalized counts using log1p for the heatmap
        self.dds.layers['log1p'] = np.log1p(self.dds.layers['normed_counts'])

        # Select significant genes for heatmap
        dds_sigs = self.dds[:, self.sigs.index]
        grapher = pd.DataFrame(
            dds_sigs.layers['log1p'].T,
            index=dds_sigs.var_names,
            columns=dds_sigs.obs_names
        )
        grapher.to_csv('grapher.csv', index=True)

        # Create a seaborn clustermap
        clustermap = sns.clustermap(grapher, z_score=0, cmap='RdYlBu_r')
        clustermap.savefig('Heatmap.png')

    def generate_volcano_plot(self):
        # Using decoupler's plot_volcano_df on self.sigs
        dc.plot_volcano_df(self.sigs, x='log2FoldChange', y='pvalue', top=30)
        plt.savefig('volcano_plot.png')

    def run_analysis(self):
        self.load_data()
        self.run_deseq_stats()
        self.perform_pca()
        self.save_pca_metadata()
        self.save_pca_proportions()
        self.save_desq2_results_table()  # Save significant genes
        self.generate_clustermap()
        self.generate_volcano_plot()

# Example usage:
# python your_script.py <count_table_file> <metadata_file> <gene_colname> <control_cond> <treat_cond> <pvalue> <log2fc>
if __name__ == "__main__":
    count_table_file = sys.argv[1]
    metadata_file = sys.argv[2]
    genes_colname = str(sys.argv[3])
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
