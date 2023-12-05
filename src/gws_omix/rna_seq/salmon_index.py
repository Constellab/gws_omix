# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com


import os

from gws_core import (ConfigParams, ConfigSpecs, File, Folder, InputSpec,
                      InputSpecs, IntParam, OutputSpec, OutputSpecs, Task,
                      TaskInputs, TaskOutputs, task_decorator)

from ..base_env.omix_env_task import BaseOmixEnvHelper


@task_decorator("SalmonIndex")
class SalmonIndex(Task):
    """
    SalmonIndex class. Represents a process that wraps Salmon index tool. Mandatory to use Salmon pseudo-aligment RNAseq mapping tool (wraps in SalmonQuantMapping class).

    Configuration options
        * `threads`: Multi threading options: number of threads to use. [Default =  2]
    """

    input_specs: InputSpecs = InputSpecs({
        'cdna_file': InputSpec(File, human_name="FastaFile", short_description="cDNA fasta file (compressed in gz)")
        # 'genome_file': InputSpec(File, human_name="FastaFile", short_description="Genome fasta file (compressed in gz)"),
        # 'gtf_annotation': InputSpec(GTFFile, human_name="GTFfile", short_description="Genome annotation file (compressed in gz)"),
    })
    output_specs: OutputSpecs = OutputSpecs({
        'salmon_index_folder': OutputSpec(
            Folder, human_name="SalmonIndexFolder", short_description="Salmon index folder")
    })

    config_specs: ConfigSpecs = {
        "threads": IntParam(default_value=2, min_value=1, short_description="Number of threads [Default =  2] ")
    }

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        thread = params["threads"]
        # annotation = inputs["gtf_annotation"]
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

        shell_proxy = BaseOmixEnvHelper.create_proxy(self.message_dispatcher)

        shell_proxy.run(cmd)

        result_folder = Folder(os.path.join(shell_proxy.working_dir, "salmon_index"))
        return {"salmon_index_folder": result_folder}
