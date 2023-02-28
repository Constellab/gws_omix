# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (File, InputSpec, IntParam, OutputSpec, TaskInputs, TaskOutputs,
                      task_decorator)
from gws_core.config.config_types import ConfigParams, ConfigSpecs
from gws_core.io.io_spec import InputSpec, OutputSpec
from gws_core.io.io_spec_helper import InputSpecs, OutputSpecs

from ..base_env.omix_env_task import BaseOmixEnvTask
#from ..file.fastq_file import FastqFile
from ..file.salmon_index_result_folder import SalmonIndexResultFolder
#from ..file.salmon_quant_mapping_output_file import SalmonQuantMappingFile

from ..file.fastq_folder import FastqFolder

from ..file.salmon_reads_quantmerge_output_file import \
    SalmonReadsQuantmergeOutputFile
from ..file.salmon_tpm_quantmerge_output_file import \
    SalmonTpmQuantmergeOutputFile


@task_decorator("SalmonQuantMapping")
class SalmonQuantMapping(BaseOmixEnvTask):
    """
    SalmonQuantMapping class. Represents a process that wraps Salmon pseudo-aligment RNAseq mapping program (salmon quant).

    Configuration options
        * `threads`: Multi threading options: number of threads to use (min=2). [Default =  2].
    """

    input_specs: InputSpecs = {
        'fastq_folder': InputSpec(FastqFolder, human_name="FastqFolder", short_description="Fastq Folder containing all RNA-sequencing file"),
        'metadata_file': InputSpec(File, human_name="MetadataFile", short_description="File containing metadata information link to file and sample (see Gencovery hub documentation at: https://hub.gencovery.com/bricks/gws_omix/latest/doc/use-cases/undefined )"),
        # 'fastq_files_forward': InputSpec(FastqFile, human_name="", short_description=""),
        # 'fastq_files_reverse': InputSpec(FastqFile, human_name="", short_description=""),
        'salmon_genome_index': InputSpec(SalmonIndexResultFolder, human_name="SalmonIndex", short_description="Folder containing Salmon index files"),
    }
    output_specs: OutputSpecs = {
        'salmon_quant_tpm':
        OutputSpec(
            SalmonTpmQuantmergeOutputFile, human_name="TpmCountFile",
            short_description="Tpm genes expression counting files merged in one (see salmon quantmerge documentation)"),
        'salmon_quant_raw':
        OutputSpec(
            SalmonReadsQuantmergeOutputFile, human_name="RawCountFile",
            short_description="Raw count genes expression counting files merged in one (see salmon quantmerge documentation)")
    }
    config_specs: ConfigSpecs = {
        'threads': IntParam(default_value=2, min_value=1, short_description="Number of threads [Default =  2] ")
    }

    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        #path = self._get_output_file_path(inputs)
        result_file_tpm = SalmonTpmQuantmergeOutputFile()
        result_file_raw = SalmonReadsQuantmergeOutputFile()
        result_file_tpm.path = os.path.join(self.working_dir, "salmon_quantmerge.tpm_count.txt")
        result_file_raw.path = os.path.join(self.working_dir, "salmon_quantmerge.raw_count.txt")
        return {
            "salmon_quant_tpm": result_file_tpm,
            "salmon_quant_raw": result_file_raw
        }

    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        thread = params["threads"]
        fq_folder = inputs["fastq_folder"]
        metadata = inputs["metadata_file"]
        #reads_forward = inputs["fastq_files_forward"]
        #reads_reverse = inputs["fastq_files_reverse"]
        transcripts_index = inputs["salmon_genome_index"]

        script_file_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [
            " bash ",
            os.path.join(script_file_dir, "./sh/salmon_quant_mapping.sh"),
            fq_folder.path,
            thread,
            metadata.path,
            transcripts_index.path
        ]

        return cmd

    # def _get_fastq_file_name(self, inputs):
    #     fastq = inputs["fastq_files_forward"]
    #     return os.path.basename(fastq.path) + ".Salmon_counting.csv"

    # def _get_output_file_path(self, inputs):
    #     return os.path.join(
    #         self.working_dir,
    #         self._get_fastq_file_name(inputs)
    #     )
