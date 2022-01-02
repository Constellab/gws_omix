# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com


import os

from gws_core import (ConfigParams, File, StrParam, TaskInputs, TaskOutputs,
                      task_decorator)

from ..base_env.omix_env_task import BaseOmixEnvTask
from ..file.salmon_reads_quantmerge_output_file import \
    SalmonReadsQuantmergeOutputFile
from ..file.salmon_tpm_quantmerge_output_file import \
    SalmonTpmQuantmergeOutputFile


@task_decorator("SalmonQuantMerge")
class SalmonQuantMerge(BaseOmixEnvTask):
    """
    Salmon quantmerge class. Represents a process that wraps Salmon quantmerge program.

    Allowed to obtain tabular files which includes a merging of previous salmon quant output files (raw reads count or normalized TPM)

    # sample.raw_reads_output.tsv
    Name    sample_1    sample_2    sample_3
    gene_1  <Raw reads count>   <Raw reads count>   <Raw reads count>
    gene_2  <Raw reads count>   <Raw reads count>   <Raw reads count>
    gene_3  <Raw reads count>   <Raw reads count>   <Raw reads count>

    # sample.tpm_output.tsv
    Name    sample_1    sample_2    sample_3
    gene_1  <Raw reads count>   <Raw reads count>   <Raw reads count>
    gene_2  <Raw reads count>   <Raw reads count>   <Raw reads count>
    gene_3  <Raw reads count>   <Raw reads count>   <Raw reads count>

    """

    output_specs = {
        'salmon_quant_tpm_file': (SalmonTpmQuantmergeOutputFile,),
        'salmon_quant_raw_reads_file': (SalmonReadsQuantmergeOutputFile,)
    }
    config_specs = {
        "experiment_name":
            StrParam(default_value="Current_experiment", short_description=" Output file names "),
    }

    def gather_output(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        return {
            "salmon_quant_tpm_file": File(path=self.output_tpm_file),
            "salmon_quant_raw_reads_file": File(path=self.output_raw_count_file)
        }

    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        brick_dir = self.get_brick_dir("omix")
        bin_file = os.path.join(brick_dir, "bin", "Salmon")

        exp_name = self.get_param("experiment_name")
        self.output_tpm_file = os.path.join(exp_name, ".tpm_output.tsv")
        self.output_raw_count_file = os.path.join(exp_name, ".raw_reads_output.tsv")
        cmd = [
            bin_file,
            " quantmerge --quants *.RNAseq_mapping.genome.Salmon_counting --column tpm -o ", self.output_tpm_file,
            " ; ",
            bin_file,
            " quantmerge --quants *.RNAseq_mapping.genome.Salmon_counting --column numreads -o ", self.output_raw_count_file,
            " ; "
        ]

        return cmd
