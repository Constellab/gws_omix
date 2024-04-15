

import os

from gws_core import (ConfigParams, ConfigSpecs, File, Folder, InputSpec,
                      InputSpecs, IntParam, OutputSpec, OutputSpecs,
                      ShellProxy, Table, Task, TaskInputs, TaskOutputs,
                      task_decorator)

from ..base_env.htseq_env_task import HTSeqShellProxyHelper


@task_decorator("HTSeqCount",  human_name="HTSeqCount",
                short_description="Genes expression count based on bam files containing folder", hide=True)
class HTSeqCount(Task):
    """
    HTSeqCount class.

    This task generate genes expression count merged files with bam containing folder. This output is compatible with DESEQ2 task.

    """
    input_specs: InputSpecs = InputSpecs({
        'bam_folder': InputSpec(Folder),
        'gtf_file': InputSpec(File)
    })
    output_specs: OutputSpecs = OutputSpecs({
        'gene_expression_file': OutputSpec(File),
        'unmapped_stats': OutputSpec(Table)
    })
    config_specs: ConfigSpecs = {
        "threads": IntParam(default_value=2, min_value=2, short_description="Number of threads")
    }

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        bam_folder: Folder = inputs["bam_folder"]
        annotation_file: File = inputs["gtf_file"]
        thrd = params["threads"]
        script_file_dir = os.path.dirname(os.path.realpath(__file__))

        shell_proxy = HTSeqShellProxyHelper.create_proxy(self.message_dispatcher)

        outputs = self.run_cmd_lines(shell_proxy,
                                     script_file_dir,
                                     bam_folder.path,
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
                      ) -> str:

        cmd_1 = [
            " . ",
            os.path.join(script_file_dir, "./sh/1_htseq_count.sh"),
            bam_folder_path,
            annotation_file,
            thrd
        ]
        self.log_info_message("HTSeq gene expression counting")
        res = shell_proxy.run(cmd_1)
        if res != 0:
            raise Exception("First step did not finished")
        self.update_progress_value(90, "Done")

        # This script perform Qiime2 demux , quality assessment
        cmd_2 = [
            ".",
            os.path.join(script_file_dir, "./sh/2_htseq_count_output_formating.sh")
        ]
        self.log_info_message("Formating output files")
        res = shell_proxy.run(cmd_2)
        if res != 0:
            raise Exception("Second step did not finished")
        self.update_progress_value(100, "Done")

        return shell_proxy.working_dir

    def outputs_annotation(self, output_folder_path: str) -> TaskOutputs:
        result_file = File(os.path.join(output_folder_path, "merged.htseq-count.txt"))

        # result_file_stat = Table()
        # result_file_stat.path = os.path.join(output_folder_path, "unmapped_stats.txt")
        return {
            "gene_expression_file": result_file,
            "unmapped_stats": None
        }
