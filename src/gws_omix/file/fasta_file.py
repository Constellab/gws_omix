# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from gws_core import (ConfigParams, File, ListParam, TextView,
                      resource_decorator, view)

from ..base_env.omix_env_task import create_omix_conda_env


@resource_decorator("FastaFile", human_name="FastaFile", short_description="Fasta File")
class FastaFile(File):
    ''' Fasta file class'''

    @view(human_name="Head And Tail TextView", view_type=TextView,
          short_description="View of the DeepEC output file first and last lines as raw text")
    def view_head_as_raw_text(self, params: ConfigParams) -> dict:
        cmd = ["head ", self.path, " ; tail ", self.path]
        shell_proxy = create_omix_conda_env()
        text = shell_proxy.check_output(cmd)
        return TextView(text)

    @view(human_name="Get Fasta Sequence",
          view_type=TextView,
          short_description="Gives a gene fasta sequence for queried genes",
          specs={"genes": ListParam(default_value=[])}
          )
    def view_query_gene_hits_as_text(self, params: ConfigParams) -> dict:  # max values = 100 gene ids
        tab = []
        genes = params["genes"]
        for gene in genes:
            cmd = ["seqtk subseq ", self.path, " <( echo ", gene, " | sed 's/ //g' )  "]
            shell_proxy = create_omix_conda_env()
            text = shell_proxy.check_output(cmd)
            tab.append(text)
        text = "\n".join(tab)
        return TextView(text)
