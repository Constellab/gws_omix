import os
from typing import Tuple
import numpy as np          # ← ajout
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scipy.cluster.hierarchy as sch
from gws_core import (CondaShellProxy, ConfigParams, ConfigSpecs, File,
                      FloatParam, InputSpec, InputSpecs, IntParam, OutputSpec,
                      OutputSpecs, PlotlyResource, ShellProxy, StrParam, Table,
                      TableImporter, Task, TaskFileDownloader, TaskInputs,
                      TaskOutputs, task_decorator)
from gws_omix.base_env.pydesq2_env_task import Pydesq2ShellProxyHelper


@task_decorator("pyDESeq2DifferentialAnalysis", human_name="pyDESeq2 pairwise differential analysis",
                short_description="Compute differential analysis using pyDESeq2 python package (pairwise comparison)")
class Pydesq2(Task):
    """
    - PyDESeq2, a Python implementation of the DESeq2 method originally developed in R, is a versatile tool for conducting differential expression analysis (DEA) on bulk RNA-seq data.<br>
    This reimplementation provides similar—but not identical—results: it achieves higher model likelihood and enables faster performance on large datasets.<br>

    - Normalization:<br>
    &nbsp;&nbsp;&nbsp;&nbsp;• <strong>Before</strong> differential expression analysis: <strong>Median-of-ratios</strong> normalization (<code>dds.deseq2()</code>) is used to correct for differences in sequencing depth across samples.<br>
    &nbsp;&nbsp;&nbsp;&nbsp;• <strong>After</strong> differential expression analysis: <strong>VST</strong> (Variance Stabilizing Transformation) is applied for exploratory analyses such as PCA, clustering, and heatmaps.<br>

    - The Wald test is used to compare two conditions.<br>
    The null hypothesis of the Wald test states that for each gene, there is no differential expression between two sample groups (e.g., treated vs. control).<br>
    If the p-value is small (e.g., p < 0.05), the null hypothesis is rejected, suggesting there is only a 5% chance that the observed difference occurred by random chance.<br>
    However, when testing many genes, a number of non-differentially expressed genes may still appear significant due to random chance (false positives).<br>

    - Results are filtered by p-value and |log2FoldChange|.<br>


    """

    input_specs: InputSpecs = InputSpecs({
        'count_table_file': InputSpec(File, human_name="count table matrix", short_description="count table matrix"),
        'metadata_file': InputSpec(File, human_name="metadata file ", short_description="tsv metadata file"),

    })

    output_specs: OutputSpecs = OutputSpecs({
        'table_1': OutputSpec(Table, human_name="Differantial expression result", short_description="Differential expression results ,providing a summary of genes showing significant changes in expression levels between two conditions."),
        'plotly_result': OutputSpec(
            PlotlyResource, human_name="Principal Component Analysis",
            short_description="Show the difference after dimensional reduction via principal component analysis."),

        'plotly_result_heatmap': OutputSpec(
            PlotlyResource, human_name="Interactive Average Hierarchical Clustering Heatmap",
            short_description=" The average linkage method and the Euclidean distance metric was used for hierarchical clustering. Utilizing scipy's hierarchical clustering, the code groups genes and samples, and then rearranges the original DataFrame based on the clustering outcomes. The resulting heatmap, created using Plotly, visually represents the reordered data, making it easier to discern patterns and relationships within the gene expression dataset."),

        'plotly_result_volcanoplot': OutputSpec(
            PlotlyResource, human_name="Interactive Volcano plot",
            short_description="This plot permit to visualize the relationship between the log2 fold change and adjusted p-values for each gene. The color scale represents log2 fold change values, and the size of the points is controlled for better visibility. The resulting plot, titled 'Volcano Plot,' provides insights into gene expression changes and their statistical significance."),

        'file_1': OutputSpec(File, human_name="heatmap", short_description="displaying average expression levels across different groups with rows representing individual genes ensembl id and columns representing samples. The hierarchical clustering dendrograms are typically displayed on the side of the heatmap, showing the relationships between samples based on their similarity in expression profiles."),
    })

    config_specs: ConfigSpecs = ConfigSpecs({
        "genes_colname": StrParam(short_description="Column name containing gene ids in expression matrix "),
        "control_condition": StrParam(short_description="normal_condition"),
        "unnormal_condition": StrParam(short_description="unnormal_condition"),
        "pvalue_value": FloatParam(default_value=0.05, min_value=0.05, short_description="pvalue_value"),
        "log2FoldChange_value": FloatParam(default_value=0.5, min_value=0.5, short_description="log2FoldChange value"),
    })

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_pydesq2_count.py"
    )

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """
        count_table_file: File = inputs['count_table_file']
        metadata_file: File = inputs['metadata_file']
        genes_colname = params["genes_colname"]
        control_condition = params["control_condition"]
        unnormal_condition = params["unnormal_condition"]
        pvalue_value = params["pvalue_value"]             # padj threshold
        log2FoldChange_value = params["log2FoldChange_value"]

        shell_proxy: ShellProxy = Pydesq2ShellProxyHelper.create_proxy(self.message_dispatcher)

        cmd = (
            f"python3 {self.python_file_path} {count_table_file.path} {metadata_file.path} "
            f"{genes_colname} {control_condition} {unnormal_condition} {pvalue_value} {log2FoldChange_value}"
        )
        result = shell_proxy.run(cmd, shell_mode=True)

        if result != 0:
            raise Exception("An error occured during the execution of the script.")

        # fichiers de sortie du script python
        heatmap_file_name      = os.path.join(shell_proxy.working_dir, "Heatmap.png")
        pydesq2_results_file   = os.path.join(shell_proxy.working_dir, "pydesq2_results_table.csv")
        pca_prop_file          = os.path.join(shell_proxy.working_dir, "pca_proportions.csv")
        pca_meta_file          = os.path.join(shell_proxy.working_dir, "pca_metadata.csv")
        grapher_file_path      = os.path.join(shell_proxy.working_dir, "grapher.csv")

        pydesq2_results_table = TableImporter.call(
            File(pydesq2_results_file),
            {'delimiter': ',', 'header': 0, 'file_format': 'csv', 'index_column': 0}
        )

        plolty_resource        = self.build_plotly(pca_meta_file, pca_prop_file)
        plolty_resource_heatmap = self.build_plotly_heatmap(grapher_file_path)
        plolty_resource_volcanoplot = self.build_plotly_volcanoplot(
            pydesq2_results_file,
            genes_colname,
            padj_thr=pvalue_value,
            log2fc_thr=log2FoldChange_value,
        )

        return {
            'file_1': File(heatmap_file_name),
            'plotly_result': plolty_resource,
            'table_1': pydesq2_results_table,
            'plotly_result_heatmap': plolty_resource_heatmap,
            'plotly_result_volcanoplot': plolty_resource_volcanoplot,
        }


# PCA function
    def build_plotly(self, pca_file_path, pca_proportion_file_path) -> PlotlyResource:
        # Read saved files
        data = pd.read_csv(pca_file_path)
        proportion_data = pd.read_csv(pca_proportion_file_path, index_col=False)

        # Create the scatter plot
        fig = px.scatter(
            data,
            x='PC1',
            y='PC2',
            color='Condition',
            text='Sample',
            size_max=10,
            opacity=0.7,
            title='PCA Plot',
            labels={'PC1': f'PC1 ({proportion_data["PC1 Proportion"].iloc[0] * 100:.2f}%)',
                    'PC2': f'PC2 ({proportion_data["PC2 Proportion"].iloc[0] * 100:.2f}%)'},
        )

        # Show the plot
        a = PlotlyResource(fig)
        a.name = "Interactive PCA"
        return a


# heatmap function

    def build_plotly_heatmap(self, grapher_file_path) -> PlotlyResource:
        # Read saved files
        # Specify the first column as the index
        data = pd.read_csv(grapher_file_path, index_col=0)

        # Extract numeric data for clustering
        # Assuming the first column is non-numeric (e.g., gene IDs)
        numeric_data = data.drop(data.columns[0], axis=1)

        # Perform hierarchical clustering
        row_linkage = sch.linkage(numeric_data.values,
                                  method='average', metric='euclidean')
        col_linkage = sch.linkage(numeric_data.values.T,
                                  method='average', metric='euclidean')

        row_order = sch.dendrogram(row_linkage, no_plot=True)['leaves']
        col_order = sch.dendrogram(col_linkage, no_plot=True)['leaves']

        # Reorder the DataFrame
        data_reordered = data.iloc[row_order, col_order]

        # Create a heatmap
        heatmap = go.Figure(data=go.Heatmap(
            z=data_reordered.values,
            x=data_reordered.columns,
            y=data_reordered.index,
            colorscale='RdYlBu'
        ))

        # Set the size of the plot
        heatmap.update_layout(width=800, height=800)

        b = PlotlyResource(heatmap)
        b.name = "Interactive Average Hierarchical Clustering Heatmap"
        return b

# Volcano plot function
    # Volcano plot function
    def build_plotly_volcanoplot(
        self,
        pydesq2_results_file_name: str,
        genes_colname: str,
        *,
        padj_thr: float,
        log2fc_thr: float,
    ) -> PlotlyResource:
        """Interactive volcano plot: -log10(padj) vs log2FC, seuils identiques à l’analyse."""

        sigs = pd.read_csv(pydesq2_results_file_name)

        # -log10(padj)
        min_nonzero = sigs.loc[sigs["padj"] > 0, "padj"].min()
        sigs["padj_plot"] = sigs["padj"].replace(0, min_nonzero)
        sigs["padj_plot"] = -np.log10(sigs["padj_plot"])

        # colonne booléenne significativité
        sigs["signif"] = (sigs["padj"] < padj_thr) & (sigs["log2FoldChange"].abs() > log2fc_thr)

        fig = px.scatter(
            sigs,
            x="log2FoldChange",
            y="padj_plot",
            color="signif",
            hover_data=[genes_colname, "padj", "log2FoldChange"],
            opacity=0.7,
            size_max=10,
            color_discrete_map={True: "#d62728", False: "#636efa"},
            title="Volcano plot : -log10(padj) vs. log2FC",
            labels={
                "log2FoldChange": "Log₂ Fold Change",
                "padj_plot": "-log₁₀(adjusted p-value)",
                "signif": f"padj<{padj_thr} & |log2FC|>{log2fc_thr}",
            },
        )

        # lignes de seuil
        fig.add_hline(y=-np.log10(padj_thr), line_dash="dot", line_color="grey",
                      annotation_text=f"padj = {padj_thr}", annotation_position="top left")
        fig.add_vline(x=log2fc_thr,  line_dash="dot", line_color="grey")
        fig.add_vline(x=-log2fc_thr, line_dash="dot", line_color="grey")

        c = PlotlyResource(fig)
        c.name = "Interactive Volcano Plot"
        return c
