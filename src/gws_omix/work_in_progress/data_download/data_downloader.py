# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

#import csv
#import json
import os

from gws_core import (File, InputSpec, IntParam, MetadataTable,
                      MetadataTableImporter, OutputSpec, Settings, StrParam,
                      Table, TableImporter, TableRowAnnotatorHelper,
                      TaskInputs, TaskOutputs, task_decorator)
from gws_core.config.config_types import ConfigParams, ConfigSpecs
from gws_core.io.io_spec import InputSpec, OutputSpec
from gws_core.io.io_spec_helper import InputSpecs, OutputSpecs
from gws_core.resource.resource_set import ResourceSet

from ..base_env.data_download_env import DataDownloadEnvTask
from ..file.fastq_folder import FastqFolder

#import re


@task_decorator("SequencingDataDownloader", human_name="Sequencing Data Downloader",
                short_description="This task uses allow to download sequencing dataset from NCBI SRA and EMBL-EBI ENA databases")
class SequencingDataDownloader(DataDownloadEnvTask):
    """
    TrimeGalore class.

    Configuration options
        "quality": IntParam(default_value=20, min_value=1, max_value=40, short_description="Trim low-quality ends from reads in addition to adapter removal. Quality score in Phred (see https://en.wikipedia.org/wiki/Phred_quality_score) [Default: 20]"),
        "paired_end": StrParam(default_value="YES",allowed_values=["YES","NO"], short_description="Sequencing strategy: paired-end file (forward reads (sample_XXX_1.fq) and reverse reads (sample_XXX_2.fq)) or not--paired. [Default: YES]"),
        "threads": IntParam(default_value=2, min_value=1, short_description="Number of threads"),
        "filter_reads_with_unsequenced_nucl": IntParam( default_value=100, short_description="The total number of Ns (as integer) a read may contain before it will be removed altogether. In a paired-end setting, either read exceeding this limit will result in the entire pair being removed from the trimmed output files."),
        "output_dir": StrParam(default_value="filtered_reads", short_description="Name of the output directory, which contained trimmed files and report files")

    """
    input_specs: InputSpecs = {
        'ID_list': InputSpec(FastqFolder, human_name="", short_description="")
    }
    output_specs: OutputSpecs = {
        "cleaned_fastq": OutputSpec(ResourceSet, human_name="", short_description=""),
        "cleaning_stats": OutputSpec(File, human_name="", short_description="")
    }
    config_specs: ConfigSpecs = {
        "db_type":
        StrParam(
            allowed_values=["NCBI", "EMBL"],
            short_description="Choose the database [Respectively, options : NCBI, EMBL]"),
        "sequeuncing_type":
        StrParam(
            allowed_values=["paired-end", "single-end"],
            short_description="Type of sequencing strategy [Respectively, options : paired-end, single-end]. Default = paired-end")
    }

    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        seq_type = params["sequencing_type"]
        filtered_folder = FastqFolder()
        stats_file = File()
        resource_set: ResourceSet = ResourceSet()
        resource_set.name = "Set of Fastq Folder(s)"

        if seq_type == "paired-end":
            singelton_folder = FastqFolder()
            filtered_folder.path = os.path.join(self.working_dir, "filtered_fastq_folder")
            filtered_folder.name = "Filtered fastq folder"
            singelton_folder.path = os.path.join(self.working_dir, "singleton_fastq_folder")
            singelton_folder.name = "Singelton fastq folder"
            stats_file.path = os.path.join(self.working_dir, "trimming_report.txt")

            resource_set.add_resource(filtered_folder)
            resource_set.add_resource(singelton_folder)
        else:
            filtered_folder.path = os.path.join(self.working_dir, "filtered_fastq_folder")
            filtered_folder.name = "Filtered fastq folder"
            stats_file.path = os.path.join(self.working_dir, "trimming_report.txt")

            resource_set.add_resource(filtered_folder)

        return {
            "cleaned_fastq": resource_set,
            "cleaning_stats": stats_file
        }

    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        fq_folder = inputs["fastq_folder"]
        thread = params["threads"]
        qual = params["quality"]
        seq_type = params["sequencing_type"]
        #hard_trim_5 = params["hard_trimming_5prime"]
        #hard_trim_3 = params["hard_trimming_3prime"]
        max_unsequenced_nucleotides = params["filter_reads_with_unsequenced_nucl"]
        min_sz = params["min_size"]

        #trimgalore_dir = params["output_dir"]

        #self._output_dir_path = self._get_output_dir_path(trimgalore_dir)

        script_file_dir = os.path.dirname(os.path.realpath(__file__))

        if seq_type == "paired-end":
            fwd = params["forward_file_differentiator"]
            rvs = params["reverse_file_differentiator"]
            script_file_dir = os.path.dirname(os.path.realpath(__file__))
            cmd = [
                "bash",
                os.path.join(script_file_dir, "./sh/trim_galore_cmd.sh"),
                fq_folder.path,
                fwd,
                rvs,
                thread,
                qual,
                # hard_trim_5,
                # hard_trim_3,
                max_unsequenced_nucleotides,
                min_sz
            ]
        else:
            cmd = [
                "bash",
                os.path.join(script_file_dir, "./sh/trim_galore_cmd_single.sh"),
                fq_folder.path,
                thread,
                qual,
                # hard_trim_5,
                # hard_trim_3,
                max_unsequenced_nucleotides,
                min_sz
            ]
        return cmd
