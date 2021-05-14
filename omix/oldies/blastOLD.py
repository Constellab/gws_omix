# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os

from gws.json import JSONData
from gws.typing import Path
from gws.shell import Shell
from gws.file import File
from gws.logger import Error

from .data import FastaFile, BlastFile

class Blast(Shell):
    """
    Blast class. Represents a process that wraps NCBI blast program.
    
    Configuration options
        * `taxo`: Kingdom name. Specify taxonomic groups to select a specific sub-set database (Faster) = bacteria, archaea, eukaryota, metazoa, chordata, mammalia, fungi, viridiplantae. [Default = all] = Slower.
        * `alignement_type`: Alignement type. Prot against Prot database (i.e blastp) or Translated Nucl against prot database (i.e blastx). Respectivly, options = PP, TNP. [Default: PP] ".
        * `evalue`: E-value to exclude results. Default = 0.00001 (i.e 1e-5).
        * `threads`: Multi threading options: number of threads to use (min=1, max=7). [Default =  4].
        * `idt`: Similarity/identity minimum percentage threshold to exclude results (min= 1, max= 100). [Default = 70].
        * `cov`: Coverage (see blast option -qcov_hsp_perc) minimum percentage threshold to exclude results (min= 1, max= 100). [Default = 70].

    """

    input_specs = { 'fasta_file': (File,) }
    output_specs = { 'blast_file': (File,),  'blast_file_parsed': (File,), 'blast_file_parsed_best_hit': (File,)  }
    config_specs = {
        "fasta": {"type": str, "Description": "Fasta input file" },
        "taxo": {"type": str, "default": "all", "allowed_values": ["all", "prokaryota", "eukaryota", "animals", "fungi", "plant"],  "Description": "Kingdom name : Specify kingdom to select the database (Faster) = prokaryota, eukaryota, animals, fungi or plant"},
        "alignement_type": {"type": str, "default": "PP", "Description": "Type of alignement to perform : Prot against Prot database (i.e blastp) or Translated Nucl against prot database (i.e blastx). [Respectivly, options : PP, TNP ]. Default = PP"}, #query_type -> prot -- prot, nucl -- nucl...
        "evalue": {"type": float, "default": 0.00001, "min": 0, "Description": "E-value : Default = 0.00001 (i.e 1e-5)"},
        "threads": {"type": int, "default": 4, "min": 1, "max": 7, "Description": "Number of threads"},
        "idt": {"type": int, "default": 70, "min": 1, "max": 100, "Description": "Similarity/identity minimum percentage threshold to exclude results : Default = 70"},
        "cov": {"type": int, "default": 70, "min": 1, "max": 100, "Description": "Coverage (see blast option -qcov_hsp_perc) minimum percentage threshold to exclude results : Default = 70"}
    }
    

    
    def after_command(self, stdout: str=None):
        self.output["blast_file"] = File( path=self.output_file )
        self.output["blast_file_parsed"] = File( path=self.output_file_parsed )
        self.output["blast_file_parsed_best_hit"] = File( path=self.output_file_parsed_best_hit )
        self.data["stdout"] = stdout
    
    def build_command(self) -> list:
        version = "2.11.0+"
        brick_dir = self.get_brick_dir("omix")    
        bin_dir = os.path.join(brick_dir, "bin", "blast+"+ version , "bin")
#        datab_dir = os.path.join(database_dir, "data", "blast+"+ version, "datab")
        datab_dir = os.path.join("/app/gws/user", "data", "blast+"+ version, "datab") # temporary, this need to lead to the db blastdirectory
        kingdom = self.get_param("taxo")
        fasta = self.get_param("fasta_file")
        alignement = self.get_param("alignment_type")     
        evalue = self.get_param("evalue")
        thread = self.get_param("threads")
        idt = self.get_param("idt")
        cov = self.get_param("cov")
        datab_file = os.path.join(datab_dir, kingdom + ".uniprotKB.faa")
       
        if alignement == "PP":
            bin_file = os.path.join(bin_dir, "blastp")
            self.output_file = os.path.join(self.cwd.name,"result.blastp.csv")
            self.output_file_parsed = os.path.join(self.cwd.name,"result.blastp_parsed.csv") 
            self.output_file_parsed_best_hit = os.path.join(self.cwd.name,"result.blastp_parsed.best_hit.csv")             
        else:
            self.output_file = os.path.join(self.cwd.name,"result.blastx.csv")
            self.output_file_parsed = os.path.join(self.cwd.name,"result.blastx_parsed.csv") 
            self.output_file_parsed_best_hit = os.path.join(self.cwd.name,"result.blastx_parsed.best_hit.csv")
            bin_file = os.path.join(bin_dir, "blastx")


        cmd = [ "conda activate ; ", #temporary, we need to activate conda env before using it
        bin_file,
        "-db ", datab_file,
        "-query ", fasta, 
        "-evalue ", evalue,
        "-num_threads ", thread,
        "-qcov_hsp_perc ", cov,
        "-outfmt  \"7 qaccver saccver pident qcovs qcovhsp length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore\" -show_gis -task \"blastp-fast\" -out ", self.output_file,
        " ; egrep -v \"^#\" ", self.output_file,
        "| perl -ne 'chomp; @t=split/\t/; if($t[2] >= ",idt,
        " ){ print $_,\"\\n\"; }' > ",self.output_file_parsed,
        "cat ", self.output_file_parsed,
        " | perl -ne 'chomp; @t=split/\\\t/; print ++\$h{\$t[0]},\"\\\t\",\$_,\"\\\n\"; ' | egrep \"^1\\\t\" | cut -f2- | egrep -v \"^#\" > ", self.output_file_parsed_best_hit                
        ]
               
        return cmd

