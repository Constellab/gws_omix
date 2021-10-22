# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import subprocess
from gws_core import File, resource_decorator, view, IntParam, TextView, ListParam, ShellEnvProxy

from ..base.omix_env_task import BaseOmixEnvTask

@resource_decorator("BlastECFile",
                    human_name="BlastECFile",
                    short_description="BlastEC File")
class BlastECFile(File):
    ''' BlastEC file class'''

    @view(view_type=TextView, human_name="Text View", short_description="View of the blastEC output file first and last lines as raw text")
    def view_head_as_raw_text(self, **kwargs) -> dict:
        cmd = ["head ", self.path, " ; tail ", self.path ]
        shell_proxy = ShellEnvProxy(BaseOmixEnvTask)
        text = shell_proxy.check_output(cmd)
        return TextView(data = text, **kwargs)

    @view(view_type=TextView, human_name="Get Gene Hits", short_description="Gives hits for queried genes", specs={"genes": ListParam(default_value=[])})
    def view_query_gene_hits_as_csv(self, genes=[], **kwargs) -> dict:
        tab=[]

# if [[ $(grep -L "$user2" /etc/passwd) ]]; then echo "No results for gene" gene;
#  echo "No results for gene" gene;
#else
#   cmd = ["cat ", self.path, " | grep -w ", '\\"' + gene +'\\" ' ]
#fi
        for gene in genes:
#            cmd = ["cat ", self.path, " | grep -w ", '\\"' + gene +'\\" | awk -v var = ' + gene + "\\'{}{if( var == TRUE ){ print $0} else{ print \"#: \"var\" not found\"}}\\'" ]
            cmd = ["cat ", self.path, " | grep -w ", '\\"' + gene +'\\" ' ]
            shell_proxy = ShellEnvProxy(BaseOmixEnvTask)
            text = shell_proxy.check_output(cmd) 
            # line = subprocess.check_output(
            #         text,
            #         shell=True
            #     )
            tab.append(text)
        text = "\n".join(tab)
        return TextView(data = text, **kwargs)

##


    # def view_query_gene_hits_as_csv(self, genes=[]) -> dict:

    #     tab=[]
    #     for gene in genes:
    #         cmd = ["cat ", self.path, " | grep -w ", gene ]
    #         env_cmd = BaseOmixEnvTask._format_command( cmd )
    #         line = subprocess.check_output(
    #                 env_cmd,
    #                 shell=True
    #             )
    #         tab.append(line)
    #     text = "\n".join(tab)
    #     return TextView(data = text)