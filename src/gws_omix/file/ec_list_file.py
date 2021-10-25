# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import subprocess
from gws_core import File, resource_decorator, view, IntParam, TextView, ShellEnvProxy #, CsvView

from ..base.omix_env_task import BaseOmixEnvTask

@resource_decorator("ECListFile",
                    human_name="ECListFile",
                    short_description="File of EC list")
class ECListFile(File):
    ''' EC list file class'''

    @view(human_name="Text View", view_type=TextView, short_description="View of the EC list file ")
    def view_as_raw_text(self, **kwargs) -> dict:
        cmd = ["cat ", self.path ]
        shell_proxy = ShellEnvProxy(BaseOmixEnvTask)
        text = shell_proxy.check_output(cmd)
        return TextView(data = text, **kwargs)
##

    # @view(human_name="TextView", view_type=TextView, short_description="View of the EC list file ")
    # def view_head_as_raw_text(self, **kwargs) -> dict:

    #     cmd = ["cat ", self.path, " | sort -u " ]
    #     env_cmd = BaseOmixEnvTask._format_command( cmd )

    #     csv = subprocess.check_output(
    #         env_cmd,
    #         shell=True
    #     )
    #     return TextView(data = csv)
