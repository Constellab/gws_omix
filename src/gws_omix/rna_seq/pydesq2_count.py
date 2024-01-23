import os
from gws_core import (ConfigParams, ConfigSpecs, File, ShellProxy,
                      FloatParam, InputSpec, InputSpecs, IntParam, OutputSpec, TableImporter,
                      OutputSpecs, StrParam, Task, TaskFileDownloader, Table,
                      TaskInputs, TaskOutputs, task_decorator , PlotlyResource , CondaShellProxy)

import plotly.express as px
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
        'table_1': OutputSpec(Table, human_name="differantial expression result", short_description="xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"),
        'plotly_result': OutputSpec(
            PlotlyResource, human_name="pathway_pca",
            short_description="Show the difference after dimensional reduction via principal component analysis."),
        'file_1': OutputSpec(File, human_name="heatmap", short_description="xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"),
        'file_2': OutputSpec(File, human_name="volcano plot", short_description="xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
    })

    config_specs: ConfigSpecs = {
        "genes_colname": StrParam(short_description="Column name containing gene ids in expression matrix "),
        "control_condition": StrParam(short_description="normal_condition"),
        "unnormal_condition": StrParam(short_description="unnormal_condition"),
        "padj_value": FloatParam(default_value=0.05, min_value=0.05, short_description="padj_value"),
        "log2FoldChange_value": FloatParam(default_value=0.5, min_value=0.5, short_description="FloatParam"),
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
        pydesq2_results_file_name = os.path.join(
            shell_proxy.working_dir, "pydesq2_results_table.csv")

        pydesq2_results_table = TableImporter.call(
            File(pydesq2_results_file_name),
            {'delimiter': ',', 'header': 0, 'file_format': 'csv', 'index_column': 0})


        # PCA
        pca_proportion_file_path = os.path.join(shell_proxy.working_dir, "pca_proportions.csv")
        pca_file_path = os.path.join(shell_proxy.working_dir, "pca_metadata.csv")
        plolty_resource = self.build_plotly(pca_file_path, pca_proportion_file_path)

        # return the output table
        return {
            'file_1': File(heatmap_file_name),
            'file_2': File(volcanoplot_file_name),
            'plotly_result': plolty_resource,
            'table_1': pydesq2_results_table

        }

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
        return PlotlyResource(fig)

