# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import subprocess
from gws_core import File, resource_decorator, view, IntParam, TextView, ShellEnvProxy #, CsvView

from ..base_env.omix_env_task import BaseOmixEnvTask

@resource_decorator("BAMToQuantFile",
                    human_name="BAMToQuantFile",
                    short_description="BAMToQuantFile")
class BAMToQuantFile(File):
    """BAMToQuantFile file class"""

    @view(view_type=TextView, human_name="Head And Tail TextView", short_description="View of the bam file first and last lines as raw text")
    def view_head_and_tail_as_raw_text(self, **kwargs) -> dict:
        cmd = ["samtools view", self.path, "|", "head ; ", "samtools view", self.path, "|", "tail"]
        shell_proxy = ShellEnvProxy(BaseOmixEnvTask)
        #run samtools view
        text = shell_proxy.check_output(cmd)
        return TextView(data = text, **kwargs)

    @view(view_type=TextView,human_name="Read length View", short_description="Read length distribution (samtools stats)")
    def view_read_length_as_box_plot(self) -> dict:
       #Read length: first column = x-axis, second column = y-axis (boxplot)
        cmd = ["samtools stats", self.path, "|", "grep \"^RL\" | cut -f 2- " ]
        shell_proxy = ShellEnvProxy(BaseOmixEnvTask)
        #run samtools stats
        text = shell_proxy.check_output(cmd)
        return TextView(data = text)


    @view(view_type=TextView, human_name="Nucleotides Content View", short_description="Nucleotides read content (samtools stats)")
    def view_nucl_count_as_piechart(self) -> dict:
        #Nucleotides content count (piechart)
        cmd = [" samtools stats", self.path, "|", "| egrep \"^[^#]+TC\" | cut -f 2- | ",
        "awk 'BEGIN{ print \"A\\\tC\\\tG\\\tT\\\tN\"} { for (i=1; i<=NF; ++i) sum[i] += $i; j=NF } END { for (i=1; i <= j; ++i) printf \"%s \", sum[i]; printf \"\\\n\"; }'"
        ]
        shell_proxy = ShellEnvProxy(BaseOmixEnvTask)
        #run samtools stats
        csv = shell_proxy.check_output(cmd)
        return TextView(data = csv)

##

    # @view(view_type=TextView, human_name="TextView", short_description="View of the bam file first and last lines as raw text")
    # def view_head_as_raw_text(self) -> dict:
    #     cmd = ["samtools view", self.path, "|", "head ; ", "samtools view", self.path, "|", "tail"  ]
    #     env_cmd = BaseOmixEnvTask._format_command( cmd )
    #     #run samtools view
    #     text = subprocess.check_output(
    #         env_cmd,
    #         shell=True
    #     )
    #     _view = TextView(data = text)
    #     return _view.to_dict()

