# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from io import StringIO

import pandas
from gws_core import (ConfigParams, File, ListParam, ShellProxy, Table,
                      TabularView, resource_decorator, view)

from ..base_env.omix_env_task import BaseOmixEnvTask


@resource_decorator("BlastECFile", human_name="BlastECFile", short_description="BlastEC File")
class BlastECFile(File):
    ''' BlastEC file class'''

    @view(view_type=TabularView, human_name="Head and Tail table",
          short_description="View of the the head and tail of the BlastEC file as table")
    def view_head_tail_as_table(self, params: ConfigParams) -> dict:
        if self.is_large():
            cmd = ["head ", self.path, " ; tail ", self.path]
        else:
            cmd = ["cat ", self.path]
        shell_proxy = ShellProxy(BaseOmixEnvTask)
        text = shell_proxy.check_output(cmd)
        file = StringIO(text)
        df = pandas.read_csv(file, sep="\t", dtype=str, header=None)
        t_view = TabularView()
        t_view.set_data(data=df)
        return t_view

    @view(view_type=TabularView,
          human_name="Gene hits table",
          short_description="Gives hits for queried genes",
          specs={"genes": ListParam(default_value=[])}
          )
    def view_query_gene_hits_as_csv(self, params: ConfigParams) -> dict:
        genes = params["genes"]
        tab = []
        for gene in genes:
            cmd = ["cat ", self.path, " | grep -w ", '\\"' + gene + '\\" ']
            shell_proxy = ShellProxy(BaseOmixEnvTask)
            text = shell_proxy.check_output(cmd)
            tab.append(text)
        text = "\n".join(tab)
        file = StringIO(text)
        df = pandas.read_csv(file, sep="\t")
        t_view = TabularView()
        t_view.set_data(data=df)
        return t_view
