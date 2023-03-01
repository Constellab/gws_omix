# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com


import os

from gws_core import (InputSpec, IntParam, OutputSpec, TaskInputs, TaskOutputs,
                      task_decorator, ConfigParams, ConfigSpecs, InputSpec, OutputSpec, InputSpecs, OutputSpecs)

from ..base_env.omix_env_task import BaseOmixEnvTask
from ..file.fasta_file import FastaFile
from ..file.gtf_file import GTFFile
from ..file.salmon_index_result_folder import SalmonIndexResultFolder


@task_decorator("SalmonIndex")
class SalmonIndex(BaseOmixEnvTask):
    """
    SalmonIndex class. Represents a process that wraps Salmon index tool. Mandatory to use Salmon pseudo-aligment RNAseq mapping tool (wraps in SalmonQuantMapping class).

    Configuration options
        * `threads`: Multi threading options: number of threads to use. [Default =  2]
    """

    input_specs: InputSpecs = {
        'cdna_file': InputSpec(FastaFile, human_name="FastaFile", short_description="cDNA fasta file (compressed in gz)")
        # 'genome_file': InputSpec(FastaFile, human_name="FastaFile", short_description="Genome fasta file (compressed in gz)"),
        # 'gtf_annotation': InputSpec(GTFFile, human_name="GTFfile", short_description="Genome annotation file (compressed in gz)"),
    }
    output_specs: OutputSpecs = {'salmon_index_folder': OutputSpec(
        SalmonIndexResultFolder, human_name="SalmonIndexFolder", short_description="Salmon index folder"), }
    config_specs: ConfigSpecs = {
        "threads": IntParam(default_value=2, min_value=1, short_description="Number of threads [Default =  2] ")
    }

    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        # path = self._get_output_file_path(inputs)

        # print("xxxxx")
        # print(path)
        # print(os.listdir(os.path.join(path)))
        # print("---")
        # with open(os.path.join(path, "ref_indexing.log")) as fp:
        #     print(fp.read())
        # print("---")

        result_folder = SalmonIndexResultFolder()
        result_folder.path = os.path.join(self.working_dir, "salmon_index")
        return {"salmon_index_folder": result_folder}

    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        thread = params["threads"]
        #annotation = inputs["gtf_annotation"]
        genome_fasta = inputs["protein_file"]
        fasta_name = os.path.basename(genome_fasta.path)
        script_file_dir = os.path.dirname(os.path.realpath(__file__))

        cmd = [
            "bash",
            os.path.join(script_file_dir, "./sh/salmon_index_cmd.sh"),
            genome_fasta.path,
            # annotation.path,
            thread,
            fasta_name
        ]
        return cmd

    # def _get_fasta_file_name(self, inputs):
    #     genome_fasta = inputs["protein_file"]
    #     return os.path.basename(genome_fasta.path)

    # def _get_output_file_path(self, inputs):
    #     return os.path.join(
    #         self.working_dir,
    #         self._get_fasta_file_name(inputs) + ".salmon_index"
    #     )
