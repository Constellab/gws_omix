# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com


import os

from gws_core import (ConfigParams, ConfigSpecs, File, InputSpec, InputSpecs,
                      IntParam, OutputSpec, OutputSpecs, TaskInputs,
                      TaskOutputs, task_decorator)
from gws_omix import FastqFolder

from ..base_env.omix_env_task import BaseOmixEnvTask
from ..file.bam_to_quant_folder import BAMToQuantFolder
from ..file.star_index_folder import StarIndexFolder


@task_decorator("STARMappingTranscriptomeReads")
class STARMappingTranscriptomeReads(BaseOmixEnvTask):
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
        'star_transcriptome_index': InputSpec(StarIndexFolder, human_name="STAR_index_folder", short_description="Folder containing STAR index files"),
        'metadata_file': InputSpec(File, human_name="Metadata File", short_description="Metadata File"),
    })

    output_specs: OutputSpecs = OutputSpecs({
        'star_bam_folder':
        OutputSpec(
            BAMToQuantFolder, human_name="Mapping_output_files",
            short_description="Mapping output file(s) (BAM format) contained in a Folder")  # ,
        # 'tpm_file':
        # OutputSpec(
        #     SalmonTpmQuantmergeOutputFile, human_name="TmpCountMergedFile",
        #     short_description="Tpm genes expression counting files merged in one (see salmon quantmerge documentation)"),
        # 'raw_count_file':
        # OutputSpec(
        #     SalmonReadsQuantmergeOutputFile, human_name="RawCountMergedFile",
        #     short_description="Raw count genes expression counting files merged in one (see salmon quantmerge documentation)")
    })

    config_specs: ConfigSpecs = {
        "threads": IntParam(default_value=2, min_value=2, short_description="Number of threads [Default =  2] "),
        "memory": IntParam(default_value=6000000000, min_value=2000000000, short_description="Memory (RAM in Bytes) usage")  # ,
        # "sequencing_strategy":
        # StrParam(
        #     allowed_values=["I", "O", "M"],
        #     short_description="To specify sequencing strategy (I = inward, O = outward, M = matching). See https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype."),
        # "library_protocol":
        # StrParam(
        #     allowed_values=["S", "U"],
        #     short_description="To specify whether the protocol for the read library takes into account strand information [S = stranded U = unstranded]. See https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype.")

    }

    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        result_folder = BAMToQuantFolder()
        # tpm_file = SalmonTpmQuantmergeOutputFile()
        # raw_file = SalmonReadsQuantmergeOutputFile()
        result_folder.path = os.path.join(self.working_dir, "star_transciptome_mapping")
        result_folder.name = "STAR BAM Folder"

        # tpm_file.path = os.path.join(self.working_dir, "salmon_quantmerge.tpm_count.txt")
        # tpm_file.name = "TPM counts merged file"

        # raw_file.path = os.path.join(self.working_dir, "salmon_quantmerge.raw_count.txt")
        # raw_file.name = "Raw counts merged file"
        return {
            'star_bam_folder': result_folder  # ,
            # 'tpm_file': tpm_file,
            # 'raw_count_file': raw_file
        }

    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        fq_folder = inputs["fastq_folder"]
        genome_index = inputs["star_transcriptome_index"]
        metadata = inputs["metadata_file"]
        thread = params["threads"]
        ram = params["memory"]
        # lib_strategy = params["library_protocol"]
        # mapping_opt = params["sequencing_strategy"]
        # experiment_opt = mapping_opt + "" + lib_strategy

        script_file_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [
            " bash ",
            os.path.join(script_file_dir, "./sh/star_transcriptome_mapping.sh"),
            fq_folder.path,
            thread,
            ram,
            metadata.path,
            genome_index.path  # ,
            # experiment_opt
        ]

        return cmd

    # def _get_output_file_path(self):
    #     return os.path.join(self.working_dir, "star_transciptome_mapping")
