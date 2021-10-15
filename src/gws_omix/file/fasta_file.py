# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import subprocess
from gws_core import File, resource_decorator, view, IntParam, TextView, ListParam, ShellEnvProxy

from ..base.omix_env_task import BaseOmixEnvTask

@resource_decorator("FastaFile",
                    human_name="FastaFile",
                    short_description="Fasta File")
class FastaFile(File):
    ''' Fasta file class'''

    @view(human_name="Head And Tail TextView", view_type=TextView, short_description="View of the DeepEC output file first and last lines as raw text")
    def view_head_as_raw_text(self, **kwargs) -> dict:
       #Read length: first column = x-axis, second column = y-axis (boxplot)
        cmd = ["head ", self.path, " ; tail ", self.path]
        shell_proxy = ShellEnvProxy(BaseOmixEnvTask)
        text = shell_proxy.check_output(cmd)
        return TextView(data = text, **kwargs)


    @view(human_name="Get Fasta Sequence", 
          view_type=TextView, 
          short_description="Gives a gene fasta sequence for queried genes",
          specs={
            "genes": ListParam(default_value=[])
         })
    def view_query_gene_hits_as_csv(self, genes=[], **kwargs) -> dict: #max values = 100 gene ids
        tab=[]
        for gene in genes:
            #get fasta sequence
            cmd = ["seqtk subseq ", self.path, " <( echo ",gene," | sed 's/ //g' )  " ]
            shell_proxy = ShellEnvProxy(BaseOmixEnvTask)
            text = shell_proxy.check_output(cmd)
            line = subprocess.check_output(
                    text,
                    shell=True
                )
            tab.append(line)
        text = "\n".join(tab)
        return TextView(data = text, **kwargs)


#
    # def view_query_gene_hits_as_csv(self, genes=[]) -> dict: #max values = 100 gene ids
    #     tab=[]
    #     for gene in genes:
    #         #get fasta sequence
    #         cmd = ["seqtk subseq ", self.path, " <( echo ",gene," | sed 's/ //g' )  " ]
    #         env_cmd = BaseOmixEnvTask._format_command( cmd )
    #         line = subprocess.check_output(
    #                 env_cmd,
    #                 shell=True
    #             )
    #         tab.append(line)
    #     text = "\n".join(tab)
    #     return TextView(data = text)
