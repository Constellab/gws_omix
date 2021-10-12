# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import subprocess
from gws_core import File, resource_decorator, view, IntParam, TextView, ShellEnvProxy, ListParam

from ..base.omix_env_task import BaseOmixEnvTask

@resource_decorator("DeepECFile",
                    human_name="DeepECFile",
                    short_description="DeepEc File")
class DeepECFile(File):
    ''' DeepEc file class'''

    @view(view_type=TextView, human_name="Text View", short_description="View of the DeepEC output file first and last lines as raw text")
    def view_as_raw_text(self, **kwargs) -> dict:
        cmd = ["head ", self.path, " ; tail ", self.path ]
        shell_proxy = ShellEnvProxy(BaseOmixEnvTask)
        text = shell_proxy.check_output(cmd)
        return TextView(data = text, **kwargs)

    @view(view_type=TextView, human_name="Get Gene Hits", short_description="Gives EC numbers for queried genes", specs={"genes": ListParam(default_value=[])})
    def view_query_gene_hits_as_csv(self, genes=[], **kwargs) -> dict:
        tab=[]
        for gene in genes:
            cmd = ["cat ", self.path, " | grep -w ", gene ]
            shell_proxy = ShellEnvProxy(BaseOmixEnvTask)
            text = shell_proxy.check_output(cmd)
            line = subprocess.check_output(
                    text,
                    shell=True
                )
            tab.append(line)
        text = "\n".join(tab)
        return TextView(data = text, **kwargs)


##

    # @view(view_type=TextView, human_name="TextView", short_description="View of the DeepEC output file first and last lines as raw text")
    # def view_head_as_raw_text(self) -> dict:
    #    #Read length: first column = x-axis, second column = y-axis (boxplot)
    #     cmd_read_length = ["cat ", self.path, "|", " head ; cat ", self.path, "|", " tail ; " ]
    #     env_cmd_read_length = BaseOmixEnvTask._format_command( cmd_read_length )

    #     csv = subprocess.check_output(
    #         env_cmd_read_length,
    #         shell=True
    #     )
    #     _view = TextView(data = [csv])
    #     return _view.to_dict()


    # @view(view_type=TextView, human_name="GetGeneHits", short_description="Gives gene EC for queried genes")
    # def view_query_gene_hits_as_csv(self) -> dict:

    #     csv={}
    #     for i in genes:
    #         cmd = ["cat ", self.path, " | grep -w ", genes[i] ]
    #         env_cmd = BaseOmixEnvTask._format_command( cmd )

    #         csv[i] = subprocess.check_output(
    #                 env_cmd,
    #                 shell=True
    #             )
    #     _view = TextView(data = [csv])
    #     return _view.to_dict()
