import os
from gws_core import (ConfigParams, ConfigSpecs, File, ShellProxy,
                      FloatParam, InputSpec, InputSpecs, IntParam, OutputSpec, TableImporter,
                      OutputSpecs, StrParam, Task, TaskFileDownloader, Table,
                      TaskInputs, TaskOutputs, task_decorator , PlotlyResource , CondaShellProxy)

import scipy.cluster.hierarchy as sch
import plotly.graph_objects as go
import plotly.express as px
from typing import Tuple
import pandas as pd

from gws_omix.base_env.pydesq2_env_task import Pydesq2ShellProxyHelper


@task_decorator("pyDESeq2DifferentialAnalysis", human_name="pyDESeq2 pairwise differential analysis",
                short_description="Compute differential analysis using pyDESeq2 python package (pairwise comparison)" , hide=True)

class Pydesq2(Task):
    """
    PyDESeq2, a Python implementation of the DESeq2 method originally developed in R (click here) , is a versatile tool for conducting differential expression analysis (DEA) with bulk RNA-seq data.
    This re-implementation yields similar, but not identical, results: it achieves higher model likelihood, allows speed improvements on large datasets.
    By implementing Wald tests, PyDESeq2 enables users to statistically evaluate the significance of these expression differences, providing a robust framework for unraveling the nuanced relationships between genes  in RNA-seq studies.

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
        'file_2': OutputSpec(File, human_name="volcano plot", short_description="Top 30 statistical significance ( represented as p-values) on the y-axis against fold change values on the x-axis for each feature (genes) in a dataset")
    })

    config_specs: ConfigSpecs = {
        "genes_colname": StrParam(short_description="Column name containing gene ids in expression matrix "),
        "control_condition": StrParam(short_description="normal_condition"),
        "unnormal_condition": StrParam(short_description="unnormal_condition"),
        "padj_value": FloatParam(default_value=0.05, min_value=0.05, short_description="padj_value"),
        "log2FoldChange_value": FloatParam(default_value=0.5, min_value=0.5, short_description="log2FoldChange value"),
    }

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_pydesq2_count.py"
    )

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """
        # retrive the input table
        count_table_file: File = inputs['count_table_file']
        metadata_file: File = inputs['metadata_file']
        genes_colname = params["genes_colname"]
        control_condition = params["control_condition"]
        unnormal_condition = params["unnormal_condition"]
        padj_value = params["padj_value"]
        log2FoldChange_value = params["log2FoldChange_value"]

        # retrieve the factor param value
        shell_proxy: ShellProxy = Pydesq2ShellProxyHelper.create_proxy(self.message_dispatcher)


        # call python file
        cmd = f"python3 {self.python_file_path} {count_table_file.path} {metadata_file.path} {genes_colname} {control_condition} {unnormal_condition} {padj_value} {log2FoldChange_value}"
        result = shell_proxy.run(cmd, shell_mode=True)

        if result != 0:
            raise Exception("An error occured during the execution of the script.")

        heatmap_file_name = os.path.join(shell_proxy.working_dir, "Heatmap.png")
        volcanoplot_file_name = os.path.join(shell_proxy.working_dir, "volcano_plot.png")

        pydesq2_results_file_name = os.path.join(shell_proxy.working_dir, "pydesq2_results_table.csv")

        pydesq2_results_table = TableImporter.call(
            File(pydesq2_results_file_name),
            {'delimiter': ',', 'header': 0, 'file_format': 'csv', 'index_column': 0})


        # PCA
        pca_proportion_file_path = os.path.join(shell_proxy.working_dir, "pca_proportions.csv")
        pca_file_path = os.path.join(shell_proxy.working_dir, "pca_metadata.csv")
        plolty_resource = self.build_plotly(pca_file_path, pca_proportion_file_path)

        #heatmap
        grapher_file_path = os.path.join(shell_proxy.working_dir, "grapher.csv")
        plolty_resource_heatmap = self.build_plotly_heatmap(grapher_file_path)

        # Volcanoplot
        plolty_resource_volcanoplot = self.build_plotly_volcanoplot(pydesq2_results_file_name)


        # return the output table
        return {
            'file_1': File(heatmap_file_name),
            'file_2': File(volcanoplot_file_name),
            'plotly_result': plolty_resource,
            'table_1': pydesq2_results_table,
            'plotly_result_heatmap' : plolty_resource_heatmap,
            'plotly_result_volcanoplot' : plolty_resource_volcanoplot


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
            labels={'PC1': f'PC1 ({proportion_data["PC1 Proportion"].iloc[0] * 100:.2f}%)', 'PC2': f'PC2 ({proportion_data["PC2 Proportion"].iloc[0] * 100:.2f}%)'},
        )

        # Show the plot
        a = PlotlyResource(fig)
        a.name = "Interactive PCA"
        return a


# heatmap function
    def build_plotly_heatmap(self, grapher_file_path) -> PlotlyResource:
        # Read saved files
        data = pd.read_csv(grapher_file_path, index_col=0)  # Specify the first column as the index

        # Extract numeric data for clustering
        numeric_data = data.drop(data.columns[0], axis=1)  # Assuming the first column is non-numeric (e.g., gene IDs)

        # Perform hierarchical clustering
        row_linkage = sch.linkage(numeric_data.values, method='average', metric='euclidean')
        col_linkage = sch.linkage(numeric_data.values.T, method='average', metric='euclidean')

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
    def build_plotly_volcanoplot(self, pydesq2_results_file_name) -> PlotlyResource:
        # Read saved files
        sigs = pd.read_csv(pydesq2_results_file_name)  # Specify the first column as the index
        # Create a volcano plot using Plotly Express
        fig = px.scatter(
            sigs,
            x='log2FoldChange',
            y='padj',
            color='log2FoldChange',  # Use 'log2FoldChange' for color scale
            hover_data=['Geneid'],  # Assuming you have a 'Geneid' column
            size_max=10,
            opacity=0.7,
            color_continuous_scale='RdBu',  # Use 'RdBu' color scale
            title='Volcano Plot',
            labels={'log2FoldChange': 'Log2 Fold Change', 'padj': 'Adjusted p-value'}
        )
        # Show the plot
        c = PlotlyResource(fig)
        c.name = "Interactive Volcano Plot"
        return c