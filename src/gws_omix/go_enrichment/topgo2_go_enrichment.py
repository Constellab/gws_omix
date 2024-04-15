

import os

from gws_core import (ConfigParams, ConfigSpecs, File, InputSpec, InputSpecs,
                      IntParam, OutputSpec, OutputSpecs, ResourceSet, StrParam,
                      TableImporter, Task, TaskInputs, TaskOutputs,
                      task_decorator)

from ..base_env.topgo2_env_task import TopGO2ShellProxyHelper

# from gws_omix import GeneList, GeneUniverse


@task_decorator("TopGO2GoTerm", human_name="TopGO2_GO_enrichment",
                short_description="Performs GO term enrichment using TopGO2 R package", hide=True)
class ToGO2GoTerm(Task):
    """
    TopGO2-GoTerm class.

    A class that wraps the TopGO2 R package for performing GO enrichment analysis on gene list data.

    TopGO2 is an open-source R pakage for performing GO enrichment analysis using the Gene Ontology (GO) database.

    Alexa A, Rahnenfuhrer J (2022). topGO: Enrichment Analysis for Gene Ontology. R package version 2.50.0.

    More information here: https://bioconductor.org/packages/release/bioc/html/topGO.html

    [Mandatory]:
        - Gene list file must contains genes, one per line.

        - The gene univers file is a tab-delimited text file format where each row contains
            as a first column gene ID followed by (<TAB>) a list of associated GO term separeted by "," (GO:XXXXXX)

    Here's are examples of a Gene Universe File and Gene List File for performing GO term enrichment analysis:

            - Gene list files :
                Gene_0102
                Gene_0708
                ...
                Gene_0909


            - Gene universe file (tab separated) :
                # Gene ontology enrichment analysis
                Gene_0001   GO:0008150, GO:0009987, GO:0016192, GO:0050896
                Gene_0002
                Gene_0003   GO:0050896
                ...
                Gene_2305   GO:001234

            OR .gmt GSEA db files (see documentation here: https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29 ) :
                # Gene ontology enrichment analysis
                GO:0000001   NA   Gene_00012  Gene_00018  Gene_00019  Gene_00011
                GO:0000002      Here is a description of the gene   Gene_00735  Gene_01257
                GO:0000003   NA   Gene_00012
                ...
                GO:0050896   NA   Gene_00018 Gene_00019
    """

    input_specs: InputSpecs = InputSpecs({
        'Gene_universe':
        InputSpec(
            File,
            short_description="A gene universe file containing all annotated genes on the genomic sequence linked to their(s) GO term (see Documentation)",
            human_name="Gene Universe file"),
        'Gene_list': InputSpec(
            File, short_description="Gene list file to assess (see Documentation)", human_name="Gene List")})

    output_specs: OutputSpecs = OutputSpecs({
        'GO_term_enrichment': OutputSpec(ResourceSet)
    })
    config_specs: ConfigSpecs = {
        'Gene_universe_format': StrParam(allowed_values=["gmt_format", "topgo_format"],
                                         short_description="Gene Universe File format (see Documentation above)"),
        'Top_results':
        IntParam(
            default_value=25, min_value=1, short_description="Top X results Fisher's Text p-value based (see TOPGO documentation)")
    }

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        gene_universe = inputs["Gene_universe"]
        gene_list = inputs["Gene_list"]
        gene_universe_format = params["Gene_universe_format"]
        top_results = params["Top_results"]

        gene_universe_path = gene_universe.path
        gene_list_file_path = gene_list.path
        script_file_dir = os.path.dirname(os.path.realpath(__file__))

        shell_proxy = TopGO2ShellProxyHelper.create_proxy(
            self.message_dispatcher)

        outputs = self.run_gsea_go(shell_proxy,
                                   script_file_dir,
                                   gene_universe_path,
                                   gene_list_file_path,
                                   gene_universe_format,
                                   top_results
                                   )
        return outputs

    def run_gsea_go(self, shell_proxy: TopGO2ShellProxyHelper,
                    script_file_dir: str,
                    gene_universe: str,
                    gene_list: str,
                    gene_universe_format: str,
                    top_results: int
                    ) -> None:
        if gene_universe_format == "gmt_format":
            cmd = [
                "bash ",
                os.path.join(script_file_dir, "./sh/gmt_to_topGO_importer.sh"),
                gene_universe
            ]
            shell_proxy.run(cmd, shell_mode=True)
            topgo_gene_universe = os.path.join(
                shell_proxy.working_dir, "TopGO2_gene_universe.csv")
            cmd_2 = [
                "Rscript --vanilla ",
                os.path.join(script_file_dir, "./R/topgo2_script.R"),
                topgo_gene_universe,
                gene_list,
                "TopGO2.GO_enrichment_results.tsv",
                top_results
            ]
            shell_proxy.run(cmd_2, shell_mode=True)
        else:
            cmd_2 = [
                "Rscript --vanilla ",
                os.path.join(script_file_dir, "./R/topgo2_script.R"),
                gene_universe,
                gene_list,
                "TopGO2.GO_enrichment_results.tsv",
                top_results
            ]
            shell_proxy.run(cmd_2, shell_mode=True)

        # Resource set
        resource_table: ResourceSet = ResourceSet()
        resource_table.name = "TopGO2 GO Term Enrichment"

        path = os.path.join(shell_proxy.working_dir,
                            "TopGO2.GO_enrichment_results.tsv")
        result_file = File(path=path)
        result_file.name = "TopGO2 GO Enrichment - All results"
        resource_table.add_resource(result_file)

        path = os.path.join(shell_proxy.working_dir,
                            "TopGO2.top_results.pvalue_0.05.csv")
        table = TableImporter.call(
            File(path=path), {'delimiter': 'tab', "index_column": 0})
        table.name = "TopGO2 GO Enrichment - Top results (Fisher Test pvalue >= 0.05)"
        resource_table.add_resource(table)

        return {
            "GO_term_enrichment": resource_table
        }
