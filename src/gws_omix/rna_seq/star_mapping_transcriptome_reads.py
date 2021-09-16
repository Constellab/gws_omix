# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, IntParam, StrParam, ConfigParams, TaskInputs, TaskOutputs, Utils
from ..base.omix_env_task import BaseOmixEnvTask
from ..file.fasta_file import FastaFile
from ..file.fastq_file import FastqFile
from ..file.gtf_file import GTFFile
from ..file.bam_to_quant_file import BAMToQuantFile


@task_decorator("STARMappingTranscriptomeReads")
class STARMappingTranscriptomeReads(BaseOmixEnvTask):
    """
    STARMappingTranscriptomeReads class. Represents a process that wraps STAR program.

    Mapping is done on the whole genome sequence but only reads overlapping gene feature are kept.
    
    --> THIS BAM OUTPUT FILE !!! IS !!! COMPATIBLE WITH SALMON QUANTIFICATION TOOL (see SalmonQuantStar class) <---
    
    Configuration options
        * `threads`: Multi threading options: number of threads to use (min=1, max=7). [Default =  4].
        * `memory`: Memory (RAM in Bytes) usage in Bytes. [Default =  8000000000 bytes (=8GB)].

    """

    input_specs = {
        'fastq_files_forward': (FastqFile,),
        'fastq_files_reverse': (FastqFile,),   
        'fasta_genome_sequence': (FastaFile,), 
        'gtf_annotation': (GTFFile,)           
    }

    output_specs = {
        'STAR_bam_file': (BAMToQuantFile,)
    }

    config_specs = {
        "threads": IntParam(default_value=12, min_value=6, description="Number of threads [Default =  12] "),
        "memory": IntParam(default_value=48000000000, min_value=8000000000, description="Memory (RAM in Bytes) usage" )
    }
   
    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        result_file = File()
        result_file.path = self._get_output_file_path(params)
        return {"STAR_bam_file": result_file}
    
    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        thread = params["threads"]
        ram = params["memory"]                   
        reads_forward  = params["fastq_files_forward"]
        reads_reverse = params["fastq_files_reverse"]
           
        cmd = [
            "STAR --genomeDir ", self.working_dir,
            " --runThreadN ", thread,
            " --limitGenomeGenerateRAM ", ram,
            " -outSAMtype BAM Unsorted --quantMode TranscriptomeSAM --quantTranscriptomeBan Singleend --readFilesIn", reads_forward, reads_reverse,
            " --outReadsUnmapped None --twopassMode Basic --readFilesCommand \"gunzip -c\" --outSAMstrandField intronMotif --outSAMunmapped Within --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --outSAMattrRGline ID:GRPundef --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 ; rm Aligned.out.bam ; mv Aligned.toTranscriptome.out.bam ", self._get_output_file_path(params)
        ] 

        return cmd

    def _get_output_file_path(self, params):
        return params["fastq_files_forward"] + ".STAR_RNAseq_genome_mapping.transcriptome_reads.unsorted.bam"