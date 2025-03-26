# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com


import os

from gws_core import (ConfigParams, ConfigSpecs, File, Folder, InputSpec,
                      InputSpecs, IntParam, OutputSpec, OutputSpecs, Task,
                      TaskInputs, TaskOutputs, task_decorator)
from gws_omix import FastqFolder

from ..base_env.omix_env_task import BaseOmixEnvHelper


@task_decorator("STARMappingTranscriptomeReads", hide=True)
class STARMappingTranscriptomeReads(Task):
    """
    STARMappingTranscriptomeReads class. Represents a process that wraps STAR program.

    Mapping is done on the whole genome sequence but only reads overlapping gene feature are kept.

    --> THIS BAM OUTPUT FOLDER !!! IS !!! COMPATIBLE WITH SALMON QUANTIFICATION TOOL (see SalmonQuantStar class) <---

    Configuration options
        * `threads`: Multi threading options: number of threads to use (min=1, max=7). [Default =  2].
        * `memory`: Memory (RAM in Bytes) usage in Bytes. [Default =  4000000000 bytes (=4GB)].

    """

    input_specs: InputSpecs = InputSpecs({
        'fastq_folder': InputSpec(FastqFolder, human_name="Folder_fastq_files", short_description="Folder containing all RNAseq sequecing files (fastq.gz)"),
        'star_transcriptome_index': InputSpec(Folder, human_name="STAR_index_folder", short_description="Folder containing STAR index files"),
        'metadata_file': InputSpec(File, human_name="Metadata File", short_description="Metadata File"),
    })

    output_specs: OutputSpecs = OutputSpecs({
        'star_bam_folder':
        OutputSpec(
            Folder, human_name="Mapping_output_files",
            short_description="Mapping output file(s) (BAM format) contained in a Folder")
    })

    config_specs: ConfigSpecs = ConfigSpecs({
        "threads": IntParam(default_value=2, min_value=2, short_description="Number of threads [Default =  2] "),
        # ,
        "memory": IntParam(default_value=6000000000, min_value=2000000000, short_description="Memory (RAM in Bytes) usage")
    })

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        fq_folder: FastqFolder = inputs["fastq_folder"]
        genome_index: Folder = inputs["star_transcriptome_index"]
        metadata: File = inputs["metadata_file"]
        thread = params["threads"]
        ram = params["memory"]
        # lib_strategy = params["library_protocol"]
        # mapping_opt = params["sequencing_strategy"]
        # experiment_opt = mapping_opt + "" + lib_strategy

        script_file_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [
            " bash ",
            os.path.join(script_file_dir,
                         "./sh/star_transcriptome_mapping.sh"),
            fq_folder.path,
            thread,
            ram,
            metadata.path,
            genome_index.path  # ,
            # experiment_opt
        ]

        shell_proxy = BaseOmixEnvHelper.create_proxy(self.message_dispatcher)
        shell_proxy.run(cmd)

        result_folder = Folder(os.path.join(
            shell_proxy.working_dir, "star_transciptome_mapping"))
        result_folder.name = "STAR BAM Folder"
        return {
            'star_bam_folder': result_folder
        }
