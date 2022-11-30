# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from typing import List, Union

from gws_core import (BarPlotView, ConfigParams, File, Folder,
                      StackedBarPlotView, StrParam, StrRField, Table,
                      TableImporter, resource_decorator, view)
from gws_core.extra import TableView


@resource_decorator("BAMToQuantFolder", human_name="BAMToQuantFolder",
                    short_description="Folder containing BAM compatible with salmon Quant")
class BAMToQuantFolder(Folder):
    """BAMToQuantFolder folder class"""

    # @view(view_type=TextView, human_name="HeadTailTextView",
    #       short_description="View of the bam file first and last lines as raw text")
    # def view_head_and_tail_as_raw_text(self, **kwargs) -> dict:
    #     cmd = ["samtools view", self.path, "|", "head ; ", "samtools view", self.path, "|", "tail"]
    #     shell_proxy = create_omix_conda_env()
    #     text = shell_proxy.check_output(cmd)
    #     return TextView(data=text, **kwargs)

    # @view(view_type=TextView, human_name="ReadLengthTextView",
    #       short_description="Read length distribution (samtools stats)")
    # def view_read_length_as_raw_text(self) -> dict:
    #     cmd = ["samtools stats", self.path, "|", "grep \"^RL\" | cut -f 2- "]
    #     shell_proxy = create_omix_conda_env()
    #     text = shell_proxy.check_output(cmd)
    #     return TextView(data=text)

    # @view(view_type=TextView, human_name="NucleotideTextView",
    #       short_description="Nucleotides read content (samtools stats)")
    # def view_nucl_count_as_raw_text(self) -> dict:
    #     cmd = [
    #         " samtools stats", self.path, "|", "| egrep \"^[^#]+TC\" | cut -f 2- | ",
    #         "awk 'BEGIN{ print \"A\\\tC\\\tG\\\tT\\\tN\"} { for (i=1; i<=NF; ++i) sum[i] += $i; j=NF } END { for (i=1; i <= j; ++i) printf \"%s \", sum[i]; printf \"\\\n\"; }'"]
    #     shell_proxy = create_omix_conda_env()
    #     csv = shell_proxy.check_output(cmd)
    #     return TextView(data=csv)