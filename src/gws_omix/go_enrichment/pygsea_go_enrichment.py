# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, File, IntParam, TableImporter, Task,
                      TaskInputs, TaskOutputs, task_decorator)
from gws_core.config.config_types import ConfigSpecs
from gws_core.io.io_spec import InputSpec, OutputSpec
from gws_core.io.io_spec_helper import InputSpecs, OutputSpecs
from gws_core.resource.resource_set import ResourceSet

from ..base_env.pygsea_pip_env import PygseaPipShellProxyHelper

# from gws_omix import GeneList, GeneUniverse


@task_decorator("GseaGoTerm", human_name="GSEA_GO_enrichment",
                short_description="Performs GO term enrichment using pyGSEA library")
class GseaGoTerm(Task):
    """
    GseaGoTerm class.

    A class that wraps the PyGSEA library for performing Gene Set Enrichment Analysis (GSEA) on gene list data.

    PyGSEA is an open-source Python library for performing gene set enrichment analysis using various databases of gene sets, including the Gene Ontology (GO) database.

    More information here: https://github.com/zqfang/gseapy

    [Mandatory]:
        - Gene list file must contains genes, one per line.

        - The Gene Matrix Transposed (GMT) file format, also known as gene univers file
            which is a tab-delimited text file format where each row represents a gene set,
            and the first column is the name of the gene set, followed by a description,
            and then a list of genes that belong to that set.

    Here's are examples of a GMT file and Gene list file for performing GO term enrichment analysis:

            - Gene list files :
                Gene_0102
                Gene_0708
                ...
                Gene_0909


            - Gene universe file (tab separated) :
                # Gene ontology enrichment analysis
                Gene_0001   NA   GO:0008150  GO:0009987  GO:0016192  GO:0050896
                Gene_0002      Here is a description of the gene   GO:0008150  GO:0009987
                Gene_0003   NA   GO:0050896
                ...
                Gene_2305   NA
    """

    input_specs: InputSpecs = {
        'Gene_universe':
        InputSpec(
            File,
            short_description="A gene universe file (.gmt format) containing all annotated genes on the genomic sequence linked to their(s) GO term (see Documentation)",
            human_name="Gene Universe file"),
        'Gene_list': InputSpec(
            File, short_description="Gene list file to assess (see Documentation)", human_name="Gene List")}
    output_specs: OutputSpecs = {
        'GO_term_enrichment': OutputSpec(ResourceSet)
    }
    config_specs: ConfigSpecs = {"Top_results_number": IntParam(
        default_value=10, min_value=1, short_description="Number of the best enriched GO term to include in the top list"),
        "Threads": IntParam(default_value=2, min_value=2, short_description="Number of threads")}

    async def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        gene_universe = inputs["Gene_universe"]
        gene_list = inputs["Gene_list"]
        top_threshold = params["Top_results_number"]
        threads = params["Threads"]

        gene_universe_path = gene_universe.path
        gene_list_file_path = gene_list.path
        script_file_dir = os.path.dirname(os.path.realpath(__file__))

        shell_proxy = PygseaPipShellProxyHelper.create_proxy(self.message_dispatcher)

        outputs = self.run_gsea_go(shell_proxy,
                                   script_file_dir,
                                   gene_universe_path,
                                   gene_list_file_path,
                                   top_threshold,
                                   threads
                                   )
        return outputs

    def run_gsea_go(self, shell_proxy: PygseaPipShellProxyHelper,
                    script_file_dir: str,
                    gene_universe: str,
                    gene_list: str,
                    top_number: int,
                    thrds: int) -> None:

        cmd = [
            "python", os.path.join(script_file_dir, "./py/_gsea_cmd.py"),
            gene_list,
            gene_universe,
            top_number,
            thrds
        ]

        shell_proxy.run(cmd, shell_mode=True)
        # shell_proxy.run_with_proxy(cmd, params=ConfigParams, inputs=TaskInputs, shell_proxy=PygseaPipShellProxyHelper)
        # Resource set
        resource_table: ResourceSet = ResourceSet()
        resource_table.name = "GO Term Enrichment Results"

        path = os.path.join(shell_proxy.working_dir, "All_results.txt")
        table = TableImporter.call(File(path=path), {'delimiter': 'tab', "index_column": 0})
        table.name = "PyGSEA GO Enrichment - All results"
        resource_table.add_resource(table)
        path = os.path.join(shell_proxy.working_dir, "Top_list.txt")
        table = TableImporter.call(File(path=path), {'delimiter': 'tab', "index_column": 0})
        table.name = "PyGSEA GO Enrichment - Top results based on the normalized enrichment score (nes)"
        resource_table.add_resource(table)

        return {
            "GO_term_enrichment": resource_table
        }
