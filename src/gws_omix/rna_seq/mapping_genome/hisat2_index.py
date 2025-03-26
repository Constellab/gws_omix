# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, Folder, IntParam, Task, InputSpecs, OutputSpecs,
                      TaskInputs, TaskOutputs, task_decorator, InputSpec, OutputSpec, ConfigSpecs, File)

from .hisat2_env import Hisat2ShellProxyHelper


@task_decorator("Hisat2_Index", human_name="Hisat2_Index",
                short_description="Build genome index for HISAT2")
class Hisat2Index(Task):
    """
    Generates an index for the reference genome using HISAT2-build.
    This index is required for alignment using HISAT2.
    """

    input_specs: InputSpecs = InputSpecs({
        'genome_file': InputSpec(File, human_name="FASTA Genome",
                                 short_description="Reference genome in FASTA format")
    })

    output_specs: OutputSpecs = OutputSpecs({
        'output': OutputSpec(Folder, human_name="HISAT2 Index",
                             short_description="Generated HISAT2 index")
    })

    config_specs: ConfigSpecs = ConfigSpecs({
        "threads": IntParam(default_value=4, min_value=1, short_description="Number of threads")
    })

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run HISAT2-build to create genome index """

        # Retrieve input parameters
        fasta_file: File = inputs['genome_file']
        threads = params["threads"]

        # Create working directory
        shell_proxy = Hisat2ShellProxyHelper.create_proxy(
            self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'hisat2_index')
        os.makedirs(result_path, exist_ok=True)

        # Define the HISAT2 index prefix
        index_prefix = os.path.join(result_path, "genome_index")

        # Run HISAT2-build (no log redirection)
        hisat2_build_cmd = f"hisat2-build -p {threads} {fasta_file.path} {index_prefix}"
        res = shell_proxy.run(hisat2_build_cmd, shell_mode=True)

        # Check if index files were created
        index_files = [f for f in os.listdir(
            result_path) if f.startswith("genome_index")]
        if res != 0 or not index_files:
            raise Exception("Error occurred while building HISAT2 index.")

        # Return the folder containing the index
        folder = Folder(result_path)
        return {'output': folder}
