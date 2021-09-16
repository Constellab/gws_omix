# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, IntParam, StrParam, FloatParam, ConfigParams, TaskInputs, TaskOutputs, Utils
from ..base.omix_env_task import BaseOmixEnvTask
from ..file.fasta_file import FastaFile
from ..file.gtf_file import GTFFile

# gmap_build -t 7  -D ./ -d INDEX.e_coli_K12.genome.fna.fasta e_coli_K12.genome.fna.fasta
#This program gmapl is designed for large genomes.
#For small genomes of less than 2^32 (4 billion) bp, please run gmap instead. --> samtools faidx -> regarder taille genome dans fai
#cat e_coli_K12.genome.fna.fasta.fai
#Chromosome	4641652	12	70	71
#samtools faidx Z_tritici.IPO323.complete_genome_and_mitochondria.fna
#cat Z_tritici.IPO323.complete_genome_and_mitochondria.fna.fai | cut -f2 | awk '{res+=$0}END{print res}'	
#--npaths 
#--min-trimmed-coverage Float
# --min-identity float
#--cross-species
# gmap -t 7 -f 2 -D ./ -d INDEX.e_coli_K12.genome.fna.fasta output.2.fa 1> output.2.fa.INDEX.e_coli_K12.genome.fna.fasta.gmapl.gff3 2> output.2.fa.INDEX.e_coli_K12.genome.fna.fasta.gmapl.log
# gffread output.2.fa.INDEX.e_coli_K12.genome.fna.fasta.gmapl.gff3 -T -o output.2.fa.INDEX.e_coli_K12.genome.fna.fasta.gmapl.gff3.gtf ; grep "path1\";" output.2.fa.INDEX.e_coli_K12.genome.fna.fasta.gmapl.gff3.gtf > output.2.fa.INDEX.e_coli_K12.genome.fna.fasta.gmapl.gff3.best_hit.gtf ; rm output.2.fa.INDEX.e_coli_K12.genome.fna.fasta.gmapl.gff3.gtf ;

@task_decorator("GmapAlign")
class GmapAlign(BaseOmixEnvTask):
    """
    Gmap alignment tool class. Represents a process that wraps Gmap aligment tool. Gmap index is Mandatory to use Gmap tools.
    
    Configuration options
        * `threads`: Multi threading options: number of threads to use. [Default =  8]
        * `min-identity"`: Do not print alignments with identity less this value (default=70.0, 0.0 would means no filtering). Note that chimeric alignments will be output regardless of this filter. 
        * `min-trimmed-coverage`: Do not print alignments with trimmed coverage less than this value (default=70.0, 0.0 would means no filtering). Note that chimeric alignments will be output regardless of this filter. 
        * `max-hit-number`: Maximum number of hits to show (default 5).  If set to 1, GMAP will not report chimeric alignments, since those imply two hits. If you want a single alignment plus chimeric alignments, then set this to be 0. [Default =  5] 
        * `cross-species`: Use a more sensitive search for canonical splicing, which helps especially for cross-species alignments and other difficult cases (genome for a far-related species/family...). [Default =  No]
    """
    input_specs = {
        'uncompressed_genome_fasta_file': (FastaFile,),
        'genome_gmap_index_name': (File,),
        'cdna_or_cds_fasta_file': (FastaFile,)          
    }
    output_specs = {
        'Gmap_gtf_file': (GTFFile,)
    }
    config_specs = {
        "threads": IntParam(default_value=8, min_value=1, description="Number of threads [Default =  8] "),
        "min-identity": FloatParam(default_value=0.7, max_value=1.0, min_value=0.0, description="Do not print alignments with identity less this value (default=0.7, 0.0 would means no filtering). Note that chimeric alignments will be output regardless of this filter. "),
        "min-trimmed-coverage": FloatParam(default_value=0.7, min_value=0.0, max_value=1.0, description="Do not print alignments with trimmed coverage less than this value (default=0.7, 0.0 would means no filtering). Note that chimeric alignments will be output regardless of this filter. "),
        "max-hit-number": IntParam(default_value=5, min_value=0, description="Maximum number of hits to show (default 5).  If set to 1, GMAP will not report chimeric alignments, since those imply two hits. If you want a single alignment plus chimeric alignments, then set this to be 0. [Default =  5] "),
        "cross-species": StrParam(default_value="No", allowed_value=["Yes","No"], description="Use a more sensitive search for canonical splicing, which helps especially for cross-species alignments and other difficult cases (genome for a far-related species/family...). [Default =  No]"),
        "alt-start-codons": StrParam(default_value="No", allowed_value=["Yes","No"], description="Also, use the alternate initiation codons (see http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). By default, without this option, only ATG is considered an initiation codon [Default =  No]"),
        "fulllength": StrParam(default_value="No", allowed_value=["Yes","No"], description="Assume full-length protein, starting with Met (ATG codon). [Default =  No]")

    }
   
    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        result_file = File()
        result_file.path = self._get_output_file_path(params)
        return {"Gmap_gtf_file": result_file} 
   
    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        thread = params["threads"]
        idt = params["min-identity"]
        cov = params["min-trimmed-coverage"]
        hit_nbr = params["max-hit-number"]

        if params["cross-species"] == "Yes":
            crs_species = " --cross-species "
        else:
            crs_species = " "
        if params["alt-start-codons"] == "Yes":
            alt_start = " --alt-start-codons "
        else:
            alt_start = " "
        if params["cross-species"] == "Yes":
            full_lght = " --fulllength "
        else:
            full_lght = " "

        genome_fasta = params["uncompressed_genome_fasta_file"]
        genome_index = params["genome_gmap_index_name"]
        genome_fai_file = genome_fasta + ".fai" 
        fa_file = params["cdna_or_cds_fasta_file"]

        cmd = [
            "$( samtools faidx ", genome_fasta, 
            "cat", genome_fai_file,
            " | cut -f2 | awk '{res+=$0}END{if(res<4000000000){print \"gmap\"}if(res>=4000000000){print \"gmapl\"}}' ) -t ", thread,
            " -f 2 --npaths ", hit_nbr,
            " --min-identity ", idt,
            " --min-trimmed-coverage ", cov,
            " ", crs_species,
            " ", alt_start,
            " ", full_lght,            
            " -D ", self.working_dir,
            "-d ", genome_index, fa_file,
            " > tmp.gff3 ; gffread tmp.gff3 -T -o ", self._get_output_file_path(params),
            " ; rm tmp.gff3 ; "
        ]           
        return cmd

    def _get_output_file_path(self, params):
        return params["cdna_or_cds_fasta_file"]+ ".alligned_on." + params["uncompressed_genome_fasta_file"] + ".gmap_alignment.gtf"