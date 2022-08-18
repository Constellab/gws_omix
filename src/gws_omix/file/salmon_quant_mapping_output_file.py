# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from io import StringIO

import pandas
from gws_core import (ConfigParams, File, IntParam, ListParam,
                      Table, TabularView, TextView, resource_decorator, view)

from ..base_env.omix_env_task import create_omix_conda_env


@resource_decorator("SalmonQuantMappingFile", human_name="SalmonQuantMappingFile",
                    short_description="SalmonQuantMapping File")
class SalmonQuantMappingFile(File):
    ''' SalmonQuantMappingFile file class'''

    @view(view_type=TabularView, human_name="Table",
          short_description="View of the SalmonQuantMappingFile output file first and last lines as raw text")
    def view_head_tail_as_table(self, params: ConfigParams) -> dict:
        if self.is_large():
            cmd = ["head ", self.path, " ; tail ", self.path]
        else:
            cmd = ["cat ", self.path]
        shell_proxy = create_omix_conda_env()
        text = shell_proxy.check_output(cmd)
        file = StringIO(text)
        df = pandas.read_csv(file, sep="\t", dtype=str, header=None)
        t_view = TabularView()
        t_view.set_data(data=df)
        return t_view

