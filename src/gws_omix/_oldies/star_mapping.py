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
from ..file.bam_file import BAMFile


@task_decorator("STARMapping")
class STARMapping(BaseOmixEnvTask):
    """
    STAR class. Represents a process that wraps STAR program.
    
    Configuration options
        * `threads`: Multi threading options: number of threads to use (min=1, max=7). [Default =  4].
        * `memory`: Memory (RAM in Bytes) usage in Bytes. [Default =  8000000000 bytes (=8GB)].
        * `alignment_type`: Mapping file output: Choose between reads mapped on the whole genome (regardless of gene position; genome option) or filtering only reads mapped on gene coordinates (transcriptome option; needed for gene counting with SALMON (TPM expresssion values): i.e. bam unsorted). [Respectivly: genome, transcriptome ]. [Default = transcriptome"].

    """

    input_specs = {
        'fastq_files_forward': (FastqFile,),
        'fastq_files_reverse': (FastqFile,),   
        'fasta_genome_sequence': (FastaFile,), 
        'gtf_annotation': (GTFFile,)           
    }

    output_specs = {
        'STAR_bam_file': (BAMFile,)
    }

    config_specs = {
        "threads": IntParam(default_value=12, min_value=6, description="Number of threads [Default =  12] "),
        "memory": IntParam(default_value=48000000000, min_value=8000000000, description="Memory (RAM in Bytes) usage" ),
        "alignment_type": StrParam(default_value="transcriptome", allowed_values=["genome", "transcriptome"], description="Mapping file output: Choose between reads mapped on the whole genome (regardless of gene position; genome option) or filtering only reads mapped on gene coordinates (transcriptome option; needed for gene counting with SALMON (TPM expresssion values): i.e. bam unsorted). [Respectivly, options : genome, transcriptome ]. Default = transcriptome"}
    }
   
    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        self.output["STAR_bam_file"] = File( path=self.output_file )   
    
    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        thread = params["threads"]
        ram = params["memory"]
        annotation = params["gtf_annotation"]  
        mapping_opt = params["alignment_type"]
        fa_genome = params["fasta_genome_sequence"]                        
        reads_forward  = params["fastq_files_forward"]
        reads_reverse = params["fastq_files_reverse"]
        transcripts_index = params["salmon_genome_index"]  

        if mapping_opt == "genome":
            self.output_file = os.path.join(self.cwd.name,".STAR_RNAseq_whole_genome_mapping.all_genome_reads.sorted.bam")
            cmd = [
                "STAR --genomeDir ",fasta_file_path,
                "--runThreadN",thread,
                "--limitGenomeGenerateRAM",ram,
                "-outSAMtype BAM Unsorted --readFilesIn",reads_forward,reads_reverse,
                "--outReadsUnmapped None --twopassMode Basic --readFilesCommand \"gunzip -c\" --outSAMstrandField intronMotif --outSAMunmapped Within --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --outSAMattrRGline ID:GRPundef --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 ; samtools sort -@ ",thread,
                " -m ", ram,
                " -o ", self.output_file,
                " Aligned.out.bam ;  rm Aligned.out.bam "
            ]      
        elif mapping_opt == "transcriptome":
            self.output_file = os.path.join(self.cwd.name,".STAR_RNAseq_genome_mapping.transcriptome_reads.unsorted.bam")
            genome_fasta = self.get_param("fasta_genome_sequence")     
            annot = self.get_param("gtf_annotation")               
            cmd = [
                bin_file,
                "--genomeDir ./ --runThreadN",thread,
                "--limitGenomeGenerateRAM",ram,
                "-outSAMtype BAM Unsorted --quantMode TranscriptomeSAM --quantTranscriptomeBan Singleend --readFilesIn",reads_forward,reads_reverse,
                "--outReadsUnmapped None --twopassMode Basic --readFilesCommand \"gunzip -c\" --outSAMstrandField intronMotif --outSAMunmapped Within --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --outSAMattrRGline ID:GRPundef --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 ; rm Aligned.out.bam ; mv Aligned.toTranscriptome.out.bam ", self.output_file
            ] 
       

        return cmd

    def _get_output_file_path(self, params):
        return params["fastq_files_forward"] + ".Salmon_RNAseq_pseudo_mapping.Salmon_counting"