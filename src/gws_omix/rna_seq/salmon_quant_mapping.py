# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, File, Folder, InputSpec,
                      InputSpecs, IntParam, OutputSpec, OutputSpecs, Task,
                      TaskInputs, TaskOutputs, task_decorator)

from ..base_env.omix_env_task import BaseOmixEnvHelper
from ..file.fastq_folder import FastqFolder


@task_decorator("SalmonQuantMapping")
class SalmonQuantMapping(Task):
    """
    SalmonQuantMapping class. Represents a process that wraps Salmon pseudo-aligment RNAseq mapping program (salmon quant).

    Configuration options
        * `threads`: Multi threading options: number of threads to use (min=2). [Default =  2].
    """

    input_specs: InputSpecs = InputSpecs({
        'fastq_folder': InputSpec(FastqFolder, human_name="FastqFolder", short_description="Fastq Folder containing all RNA-sequencing file"),
        'metadata_file': InputSpec(File, human_name="MetadataFile", short_description="File containing metadata information link to file and sample (see Gencovery hub documentation at: https://hub.gencovery.com/bricks/gws_omix/latest/doc/use-cases/undefined )"),
        # 'fastq_files_forward': InputSpec(File, human_name="", short_description=""),
        # 'fastq_files_reverse': InputSpec(File, human_name="", short_description=""),
        'salmon_genome_index': InputSpec(Folder, human_name="SalmonIndex", short_description="Folder containing Salmon index files"),
    })
    output_specs: OutputSpecs = OutputSpecs({
        'salmon_quant_tpm':
        OutputSpec(
            File, human_name="TpmCountFile",
            short_description="Tpm genes expression counting files merged in one (see salmon quantmerge documentation)"),
        'salmon_quant_raw':
        OutputSpec(
            File, human_name="RawCountFile",
            short_description="Raw count genes expression counting files merged in one (see salmon quantmerge documentation)")
    })
    config_specs: ConfigSpecs = {
        'threads': IntParam(default_value=2, min_value=1, short_description="Number of threads [Default =  2] ")
    }

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        shell_proxy = BaseOmixEnvHelper.create_proxy(self.message_dispatcher)

        thread = params["threads"]
        fq_folder: FastqFolder = inputs["fastq_folder"]
        metadata: File = inputs["metadata_file"]
        # reads_forward = inputs["fastq_files_forward"]
        # reads_reverse = inputs["fastq_files_reverse"]
        transcripts_index: Folder = inputs["salmon_genome_index"]

        script_file_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [
            " bash ",
            os.path.join(script_file_dir, "./sh/salmon_quant_mapping.sh"),
            fq_folder.path,
            thread,
            metadata.path,
            transcripts_index.path
        ]

        shell_proxy.run(cmd)

        # path = self._get_output_file_path(inputs)
        result_file_tpm = File(os.path.join(shell_proxy.working_dir, "salmon_quantmerge.tpm_count.txt"))
        result_file_raw = File(os.path.join(
            shell_proxy.working_dir, "salmon_quantmerge.raw_count.txt"))
        return {
            "salmon_quant_tpm": result_file_tpm,
            "salmon_quant_raw": result_file_raw
        }
