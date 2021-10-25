# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import subprocess
from gws_core import File, resource_decorator, view, IntParam, TextView, ListParam, ShellEnvProxy #, CsvView

from ..base.omix_env_task import BaseOmixEnvTask


@resource_decorator("GTFFile",
                    human_name="GTFFile",
                    short_description="GTF File")
class GTFFile(File):
    ''' GTF file class'''

    @view(view_type=TextView, human_name="Text View", short_description="View of the annotation file first and last lines as raw text")
    def view_head_as_raw_text(self, **kwargs) -> dict:
        cmd = ["head ", self.path, " ; tail ", self.path ]
        shell_proxy = ShellEnvProxy(BaseOmixEnvTask)
        text = shell_proxy.check_output(cmd)
        return TextView(data = text, **kwargs)

    @view(view_type=TextView, human_name="Get Annotation features", short_description="Gives annotation for queried genes", specs={"genes": ListParam(default_value=[])})
    def view_query_gene_hits_as_csv(self, genes=[], **kwargs) -> dict:
        tab=[]
        for gene in genes:
            cmd = ["cat ", self.path, " | grep -w ", gene ]
            shell_proxy = ShellEnvProxy(BaseOmixEnvTask)
            text = shell_proxy.check_output(cmd)
            # line = subprocess.check_output(
            #         text,
            #         shell=True
            #     )
            tab.append(text)
        text = "\n".join(tab)
        return TextView(data = text, **kwargs)