# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, IntParam, StrParam, FloatParam, ConfigParams, TaskInputs, TaskOutputs, Utils, Settings, Folder
from ..base.trimfq_env_task import TrimFqEnvTask
from ..file.fastq_file import FastqFile


@task_decorator("TrimGalore")
class TrimGalore(TrimFqEnvTask):
    """
    TrimeGalore class. Represents a process that wraps NCBI blast program. This version !!! ALLOWED !!! to get EC numbers for digital twins reconstruction.
    
    Configuration options
        "quality": IntParam(default_value=20, min_value=1, max_value=40, short_description="Trim low-quality ends from reads in addition to adapter removal. Quality score in Phred (see https://en.wikipedia.org/wiki/Phred_quality_score) [Default: 20]"),
        "paired_end": StrParam(default_value="YES",allowed_values=["YES","NO"], short_description="Sequencing strategy: paired-end file (forward reads (sample_XXX_1.fq) and reverse reads (sample_XXX_2.fq)) or not--paired. [Default: YES]"),
        "threads": IntParam(default_value=2, min_value=1, short_description="Number of threads"),
        "filter_reads_with_unsequenced_nucl": IntParam( default_value=100, short_description="The total number of Ns (as integer) a read may contain before it will be removed altogether. In a paired-end setting, either read exceeding this limit will result in the entire pair being removed from the trimmed output files."),
        "output_dir": StrParam(default_value="filtered_reads", short_description="Name of the output directory, which contained trimmed files and report files")

    """
   
    input_specs = {
        'fastq_file_1': (FastqFile,),
        'fastq_file_2': (FastqFile,)
    }
    output_specs = {
        'filtered_dir': (Folder,)
    }
    config_specs = {
        "quality": IntParam(default_value=20, min_value=1, max_value=40, short_description="Trim low-quality ends from reads in addition to adapter removal. Quality score in Phred (see https://en.wikipedia.org/wiki/Phred_quality_score) [Default: 20]"),
        "paired_end": StrParam(default_value="YES",allowed_values=["YES","NO"], short_description="Sequencing strategy: paired-end file (forward reads (sample_XXX_1.fq) and reverse reads (sample_XXX_2.fq)) or not--paired. [Default: YES]"),
        "threads": IntParam(default_value=2, min_value=1, short_description="Number of threads"),
        "filter_reads_with_unsequenced_nucl": IntParam( default_value=100, short_description="The total number of Ns (as integer) a read may contain before it will be removed altogether. In a paired-end setting, either read exceeding this limit will result in the entire pair being removed from the trimmed output files."),
        "output_dir": StrParam(default_value="filtered_reads", short_description="Name of the output directory, which contained trimmed files and report files")
    }


    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        result_file = Folder()
        result_file = Folder(path=self._output_file_path)
        return {"filtered_dir": result_file} 

    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        qual  = params["quality"]
        paired  = params["paired_end"]
        max_unsequenced_nucleotides = params["filter_reads_with_unsequenced_nucl"]
        thread  = params["threads"]
        trimgalore_dir = params["output_dir"]

        fq_file_1 = inputs["fastq_file_1"]
        fq_file_1_name = os.path.basename(fq_file_1.path)
        fq_file_2 = inputs["fastq_file_2"]
        fq_file_2_name = os.path.basename(fq_file_2.path)

        self._output_dir_path = self._get_output_dir_path(trimgalore_dir)

        script_file_dir = os.path.dirname(os.path.realpath(__file__))

        cmd = [  
            "bash", 
            os.path.join(script_file_dir, "./sh/trim_galore_cmd.sh"),
            fq_file_1_name,
            fq_file_2_name,
            thread,
            qual,
            max_unsequenced_nucleotides,
            paired,
            self._output_file_path
        ]                 

        return cmd

    def _get_output_dir_path(self, output_directory_name) :
        return os.path.join(
            self.working_dir, 
            output_directory_name + "_filterted" 
        )
