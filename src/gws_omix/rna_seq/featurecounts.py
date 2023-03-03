# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, InputSpec, InputSpecs,
                      IntParam, OutputSpec, OutputSpecs, ShellProxy, Task,
                      TaskInputs, TaskOutputs, task_decorator)

from ..base_env.subread_env_task import SubreadShellProxyHelper
from ..file.bam_to_quant_folder import BAMToQuantFolder
from ..file.gtf_file import GTFFile
from ..file.salmon_reads_quantmerge_output_file import \
    SalmonReadsQuantmergeOutputFile


@task_decorator("FeatureCounts",  human_name="Feature Counts",
                short_description="Genes expression count based on bam files containing folder")
class FeatureCounts(Task):
    """
    featureCounts class.

    This task generate genes expression count merged files with bam containing folder. This output is compatible with DESEQ2 task.

    """
    input_specs: InputSpecs = {
        'bam_folder': InputSpec(BAMToQuantFolder),
        'gtf_file': InputSpec(GTFFile)
    }
    output_specs: OutputSpecs = {
        'gene_expression_file': OutputSpec(SalmonReadsQuantmergeOutputFile)  # ,
        # 'summary_stats': OutputSpec(Table)
    }
    config_specs: ConfigSpecs = {
        "threads": IntParam(default_value=2, min_value=2, short_description="Number of threads")
    }

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        bam_folder = inputs["bam_folder"]
        bam_folder_path = bam_folder.path
        annotation_file: GTFFile = inputs["gtf_file"]
        thrd = params["threads"]
        script_file_dir = os.path.dirname(os.path.realpath(__file__))

        shell_proxy = SubreadShellProxyHelper.create_proxy(self.message_dispatcher)

        outputs = self.run_cmd_lines(shell_proxy,
                                     script_file_dir,
                                     bam_folder_path,
                                     annotation_file.path,
                                     thrd
                                     )
        # Output formating and annotation

        annotated_outputs = self.outputs_annotation(outputs)

        return annotated_outputs

    def run_cmd_lines(self, shell_proxy: ShellProxy,
                      script_file_dir: str,
                      bam_folder_path: str,
                      annotation_file: str,
                      thrd: int
                      ) -> None:

        cmd_1 = [
            " . ",
            os.path.join(script_file_dir, "./sh/1_featurecounts.sh"),
            bam_folder_path,
            annotation_file,
            thrd
        ]
        self.log_info_message("featureCounts gene expression counting")
        res = shell_proxy.run(cmd_1)
        if res != 0:
            raise Exception("First step did not finished")
        self.update_progress_value(90, "Done")

        cmd_2 = [
            ".",
            os.path.join(script_file_dir, "./sh/2_featurecounts_output_formating.sh")
        ]
        self.log_info_message("Formating output files")
        res = shell_proxy.run(cmd_2)
        if res != 0:
            raise Exception("Second step did not finished")
        self.update_progress_value(100, "Done")

        output_folder_path = shell_proxy.working_dir

        return output_folder_path

    def outputs_annotation(self, output_folder_path: str) -> None:
        result_file = SalmonReadsQuantmergeOutputFile()
        # result_file_stat = Table()
        result_file.path = os.path.join(output_folder_path, "featurecounts.gene_expression_matrix.txt")
        # result_file_stat.path = os.path.join(output_folder_path, "featurecounts.summary.txt")
        return {
            "gene_expression_file": result_file  # ,
            # "summary_stats": result_file_stat
        }
