# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import subprocess

from gws_core import File, IntParam, TextView, resource_decorator, view

from ..base_env.omix_env_task import BaseOmixEnvTask


@resource_decorator("salmon_reads_quantmerge_output_file",
                    human_name="salmon_reads_quantmerge_output_file",
                    short_description="salmon_reads_quantmerge_output_file")
class SalmonReadsQuantmergeOutputFile(File):
    """salmon_Reads_quantmerge_output_file class"""

    @view(view_type=TextView, human_name="TextView",
          short_description="View of the expression file first and last lines as raw text")
    def view_head_as_raw_text(self) -> dict:
        cmd_read_length = ["cat ", self.path, "|", " head ; cat ", self.path, "|", " tail ; "]
        env_cmd_read_length = BaseOmixEnvTask._format_command(cmd_read_length)
        csv = subprocess.check_output(
            env_cmd_read_length,
            shell=True
        )
        _view = TextView(data=[csv])
        return _view.to_dict()

    @view(view_type=TextView, human_name="ExpressionInTpmPerSamples",
          short_description="Summary of genes expression (in reads mapped) per samples")
    def view_tpm_per_columns_as_box_plot(self) -> dict:
        cmd = ["cat ", self.path]
        env_cmd = BaseOmixEnvTask._format_command(cmd)
        csv = subprocess.check_output(
            env_cmd,
            shell=True
        )
        _view = TextView(data=[csv])
        return _view.to_dict()
