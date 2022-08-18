# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from gws_core import (ConfigParams, File, IntParam, ListParam,
                      TextView, resource_decorator, view)

from ..base_env.omix_env_task import create_omix_conda_env


@resource_decorator("GFF3File", human_name="GFF3File", short_description="GFF3 File")
class GFF3File(File):
    ''' GFF3 file class'''

    @view(view_type=TextView, human_name="TextView",
          short_description="View of the annotation file first and last lines as raw text")
    def view_head_as_raw_text(self, params: ConfigParams) -> dict:
        cmd = ["head ", self.path, " ; tail ", self.path]
        shell_proxy = create_omix_conda_env()
        text = shell_proxy.check_output(cmd)
        return TextView(text)

    @view(view_type=TextView, human_name="AnnotationTextView",
          short_description="Gives annotation for queried genes",
          specs={"genes": ListParam(default_value=[])}
          )
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
