# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import subprocess
from gws_core import File, resource_decorator, view, IntParam, TextView, ListParam, ShellEnvProxy, CsvView

from ..base.omix_env_task import BaseOmixEnvTask


@resource_decorator("SalmonReadsQuantmergeOutputFile",
                    human_name="Salmon Read counts Quantmerge Output File",
                    short_description="Output file from Salmon Quantmerge (expression in reads mapped)")
class SalmonReadsQuantmergeOutputFile(File):
    """salmon_Reads_quantmerge_output_file class"""

    @view(view_type=CsvView, human_name="Raw text", short_description="View of the expression file first and last lines as raw text")
    def view_head_as_raw_text(self, **kwargs) -> dict:
        cmd = ["head ", self.path, " ; tail ", self.path ]
        shell_proxy = ShellEnvProxy(BaseOmixEnvTask)
        text = shell_proxy.check_output(cmd)
        return CsvView(data = text, **kwargs)

    # @view(view_type=CsvView, human_name="Expression In Read Counts Per Samples", short_description="Summary of genes expression (in reads mapped) per samples")
    # def view_reads_counts_per_columns_as_box_plot(self) -> dict:
    #    #Plotting 1 boxplot per columns (e.g. per samples) without showing the .csv (separator= tab) file
    #     cmd = ["cat ", self.path ]
    #     shell_proxy = ShellEnvProxy(BaseOmixEnvTask)
    #     text = shell_proxy.check_output(cmd)
    #     # csv = subprocess.check_output(
    #     #     text,
    #     #     shell=True
    #     # )
    #     return TextView(data = csv)

    @view(view_type=CsvView, human_name="Get Gene Expression", short_description="Gives gene expression values for queried genes", specs={"genes": ListParam(default_value=[])})
    def view_query_gene_tpm_as_csv(self, genes=[], **kwargs) -> dict:
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
        return CsvView(data = text, **kwargs)

##

    # @view(view_type=TextView, human_name="TextView", short_description="View of the expression file first and last lines as raw text")
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
