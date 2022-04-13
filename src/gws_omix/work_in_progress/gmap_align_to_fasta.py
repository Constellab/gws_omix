# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, FloatParam, IntParam, StrParam, ConfigParams, TaskInputs, TaskOutputs, Utils
from ..base_env.omix_env_task import BaseOmixEnvTask
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

@task_decorator("GmapAlignFasta",)
class GmapAlignFasta(BaseOmixEnvTask):
    """
    GmapAlignFasta class. 
    
    Represents a task that wraps Gmap aligment tool. Gmap index is mandatory to use Gmap tools.
    
    Configuration options
        * `threads`: Multi threading options: number of threads to use. [Default =  8]
        * `min-identity"`: Do not print alignments with identity less this value (default=70.0, 0.0 would means no filtering). Note that chimeric alignments will be output regardless of this filter. 
        * `min-trimmed-coverage`: Do not print alignments with trimmed coverage less than this value (default=70.0, 0.0 would means no filtering). Note that chimeric alignments will be output regardless of this filter. 
        * `max-hit-number`: Maximum number of hits to show (default 5).  If set to 1, GMAP will not report chimeric alignments, since those imply two hits. If you want a single alignment plus chimeric alignments, then set this to be 0. [Default =  5] 
        * `cross-species`: Use a more sensitive search for canonical splicing, which helps especially for cross-species alignments and other difficult cases (genome for a far-related species/family...). [Default =  No]
    """
    input_specs = {
        'uncompressed_genome_fasta_file': (FastaFile,),
        'cdna_or_cds_fasta_file': (FastaFile,)          
    }
    output_specs = {
        'gmap_fasta_file': (FastaFile,)
    }
    config_specs = {
        "threads": IntParam(default_value=8, min_value=1, short_description="Number of threads [Default =  8] "),
        "min-identity": FloatParam(default_value=0.7, max_value=1.0, min_value=0.0, short_description="Do not print alignments with identity less this value (default=0.7, 0.0 would means no filtering). Note that chimeric alignments will be output regardless of this filter. "),
        "min-trimmed-coverage": FloatParam(default_value=0.7, min_value=0.0, max_value=1.0, short_description="Do not print alignments with trimmed coverage less than this value (default=0.7, 0.0 would means no filtering). Note that chimeric alignments will be output regardless of this filter. "),
        "max-hit-number": IntParam(default_value=5, min_value=0, short_description="Maximum number of hits to show (default 5).  If set to 1, GMAP will not report chimeric alignments, since those imply two hits. If you want a single alignment plus chimeric alignments, then set this to be 0. [Default =  5] "),
        "cross-species": StrParam(default_value="No", allowed_values=["Yes","No"], short_description="Use a more sensitive search for canonical splicing, which helps especially for cross-species alignments and other difficult cases (genome for a far-related species/family...). [Default =  No]"),
        "alt-start-codons": StrParam(default_value="No", allowed_values=["Yes","No"], short_description="Also, use the alternate initiation codons (see http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). By default, without this option, only ATG is considered an initiation codon [Default =  No]"),
        "fulllength": StrParam(default_value="No", allowed_values=["Yes","No"], short_description="Assume full-length protein, starting with Met (ATG codon). [Default =  No]")

    }
   

    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        result_file = FastaFile()
        result_file.path = self._output_file_path
        return {"gmap_fasta_file": result_file} 
   
    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        thread = params["threads"]
        idt = params["min-identity"]
        cov = params["min-trimmed-coverage"]
        hit_nbr = params["max-hit-number"]
        if params["cross-species"] == "Yes":
            crs_species = "Y" #" --cross-species "
        else:
            crs_species = "N"
        if params["alt-start-codons"] == "Yes":
            alt_start = "Y" #" --alt-start-codons "
        else:
            alt_start = "N"
        if params["cross-species"] == "Yes":
            full_lgth = "Y" #" --fulllength "
        else:
            full_lgth = "N"

        genome_fasta = inputs["uncompressed_genome_fasta_file"]
        genome_fasta_file_name = os.path.basename(genome_fasta.path)
        fa_file = inputs["cdna_or_cds_fasta_file"]
        fasta_file_name = os.path.basename(fa_file.path)
        self._output_file_path = self._get_output_file_path(fasta_file_name, genome_fasta_file_name)
        script_file_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [
            "bash", 
            os.path.join(script_file_dir, "./sh/gmap_fasta_cmd.sh"),
            thread,
            idt,
            cov,
            hit_nbr,
            crs_species,
            alt_start,
            full_lgth,
            genome_fasta.path,
            fa_file.path,
            self._output_file_path            
        ]
        return cmd

    def _get_output_file_path(self, fasta, genome):
        return os.path.join(
            self.working_dir, 
            fasta + ".alligned_on." + genome + ".gmap_alignment.fasta"
        )