# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from gws_core import (ConfigParams, File, ListParam, TextView,
                      resource_decorator, view)

from ..base_env.omix_env_task import create_omix_conda_env


@resource_decorator("DeepECFile", human_name="DeepECFile", short_description="DeepEc File")
class DeepECFile(File):
    ''' DeepEc file class'''

    @view(view_type=TextView, human_name="TextView",
          short_description="View of the DeepEC output file first and last lines as raw text")
    def view_as_raw_text(self, params: ConfigParams) -> dict:
        cmd = ["head ", self.path, " ; tail ", self.path]
        shell_proxy = create_omix_conda_env()
        text = shell_proxy.check_output(cmd)
        return TextView(text)

    @view(view_type=TextView, human_name="GeneHitsTextView", short_description="Gives EC numbers for queried genes",
          specs={"genes": ListParam(default_value=[])})
    def view_query_gene_hits_as_csv(self, params: ConfigParams) -> dict:
        tab = []
        genes = params["genes"]
        for gene in genes:
            cmd = ["cat ", self.path, " | grep -w ", gene]
            shell_proxy = create_omix_conda_env()
            text = shell_proxy.check_output(cmd)
            tab.append(text)
        text = "\n".join(tab)
        return TextView(text)
