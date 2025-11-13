# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (
    ConfigParams, Folder, IntParam, StrParam, Task, InputSpecs, OutputSpecs,
    TaskInputs, TaskOutputs, task_decorator, InputSpec, OutputSpec, ConfigSpecs, File
)

from .salmon_env import SalmonShellProxyHelper

@task_decorator("Salmon_Index", human_name="Salmon_Index",
                short_description="Build Salmon transcriptome index")
class SalmonIndex(Task):
    """
    This task builds a Salmon transcriptome index from a transcript FASTA. Provide transcript_fasta (FASTA of transcripts). Configure threads (CPU cores; default 2). The task runs:
    salmon index -t <transcripts.fasta> -i <working_dir>/salmon_index -p <threads>
    and returns output pointing to <working_dir>/salmon_index/. Use this index later with Salmon quantification tasks (salmon quant).
    """
    input_specs: InputSpecs = InputSpecs({
        'transcript_fasta': InputSpec(File, human_name="Transcript FASTA",
                                      short_description="Fasta file containing transcripts")
    })

    output_specs: OutputSpecs = OutputSpecs({
        'output': OutputSpec(Folder, human_name="Salmon Index",
                             short_description="Folder containing Salmon index files")
    })

    config_specs: ConfigSpecs = ConfigSpecs({
        'threads': IntParam(default_value=2, min_value=1, short_description="Number of threads")
    })

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """
        1) Create working dir
        2) salmon index -t <transcripts.fasta> -i <index_path> -p <threads>
        3) Return index folder
        """
        salmon_proxy = SalmonShellProxyHelper.create_proxy(self.message_dispatcher)

        transcripts: File = inputs['transcript_fasta']
        threads = params["threads"]

        result_path = os.path.join(salmon_proxy.working_dir, "salmon_index")
        os.makedirs(result_path, exist_ok=True)

        index_cmd = (
            f"salmon index "
            f"-t {transcripts.path} "
            f"-i {result_path} "
            f"-p {threads}"
        )
        print("[DEBUG] Salmon Index CMD:", index_cmd)
        res = salmon_proxy.run(index_cmd, shell_mode=True)
        if res != 0:
            raise Exception("Salmon index failed.")

        # Return the folder with index files
        return {'output': Folder(result_path)}
