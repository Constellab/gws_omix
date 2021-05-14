# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
from gws.shell import CondaShell

from gws.json import JSONData
from gws.typing import Path
from gws.shell import Shell
from gws.file import File
from gws.logger import Error
from gws.settings import Settings

class BlastEC(CondaShell):
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
    
    VERSION = "2.11.0+"
    
    input_specs = { 'fasta_file': (File,) }
    output_specs = { 'raw_blast_file': (File,), 'filtered_blast_file': (File,)  }
    config_specs = {
        "taxo": {"type": str, "default": "fungi", "allowed_values": ["all", "prokaryota", "eukaryota", "animals", "fungi", "plant"],  "Description": "Kingdom name : Specify kingdom to select the database (Faster) = prokaryota, eukaryota, animals, fungi or plant"},
        "alignment_type": {"type": str, "default": "PP", "Description": "Type of alignement to perform : Prot against Prot database (i.e blastp) or Translated Nucl against prot database (i.e blastx). [Respectivly, options : PP, TNP ]. Default = PP"}, #query_type -> prot -- prot, nucl -- nucl...
        "e_value": {"type": float, "default": 0.00001, "min": 0, "Description": "E-value : Default = 0.00001 (i.e 1e-5)"},
        "threads": {"type": int, "default": 4, "min": 1, "max": 7, "Description": "Number of threads"},
        "idt": {"type": int, "default": 70, "min": 1, "max": 100, "Description": "Similarity/identity minimum percentage threshold to exclude results : Default = 70"},
        "cov": {"type": int, "default": 70, "min": 1, "max": 100, "Description": "Coverage (see blast option -qcov_hsp_perc) minimum percentage threshold to exclude results : Default = 70"}
    }
    

    
# functions

    def gather_outputs(self, stdout: str=None):
        # execute blast cmd        
        self.output["raw_blast_file"] = File( path=self.raw_output_file ) # Blast output
        
        # execute blast parsing command and ec number retrieving 
        settings = Settings.retrieve()
        data_dir = settings.get_data_dir()
        datab_dir = os.path.join(data_dir ,"prod","omix", "blast+"+ self.VERSION )
        
        taxo = self.get_param("taxo")
        idt = self.get_param("idt")
        tab_file = os.path.join(datab_dir, taxo + ".uniprotKB.tab")        
        
        self._create_filtered_output_file(tab_file, idt)
        self.output["filtered_blast_file"] = File( path=self.filtered_output_file ) # blast output with EC number and best hit info.
        self.data["stdout"] = stdout
        
    def build_command(self) -> list:
        settings = Settings.retrieve()
        data_dir = settings.get_data_dir()
        brick_dir = settings.get_cwd()  
        bin_dir = os.path.join(brick_dir, "bin", "blast+"+ self.VERSION , "bin")
        datab_dir = os.path.join(data_dir ,"prod","omix", "blast+"+ self.VERSION ) # temporary, this need to lead to the db blastdirectory 
        taxo = self.get_param("taxo")   
        alignment = self.get_param("alignment_type")     
        evalue = self.get_param("e_value")
        thread = self.get_param("threads")
        idt = self.get_param("idt")
        cov = self.get_param("cov")
        datab_file_path = os.path.join(datab_dir, taxo + ".uniprotKB.faa")

        if alignment == "PP":
            blast_binary = "blastp"
            self.raw_output_file = os.path.join(self.cwd.name,"raw_blastp_result.tsv")
            self.filtered_output_file = os.path.join(self.cwd.name,"filtered_blastp_result.tsv")
        else:
            blast_binary = "blastx"
            self.raw_output_file = os.path.join(self.cwd.name,"raw_blastx_result.tsv")
            self.filtered_output_file = os.path.join(self.cwd.name,"filtered_blastx_result.tsv")
        
        fasta_file = self.input["fasta_file"]
        bash_file_path = os.path.join(brick_dir, "bin", "blast_ec.sh")
        cmd = [ 
            "bash",
            bash_file_path,
            blast_binary,
            datab_file_path,
            fasta_file.path,
            evalue,
            thread,
            cov,
            idt,
            self.raw_output_file,
        ]
        return cmd

    def _create_filtered_output_file(self, tabular_file, id):
        gene_ec={}
        hit_parsed={}

        with open(tabular_file, 'r') as lines: # Create dict. containing genes with their corresponding EC number(s)
            li=lines.readlines()
            for i, line in enumerate(li):
                if re.match("^#",line):
                    pass
                else:
                    li_split=line.split("\t")
                    if re.match("^$",str(li_split[7])):
                        gene_ec[li_split[0]]="NA\n"                   
                    else:
                        gene_ec[li_split[0]]=li_split[7]                    

        
        with open(self.filtered_output_file, 'a') as filtered_fp:
            with open(self.raw_output_file, 'r') as raw_fp: 
                # Create dict. containing for each lines of the blast output 
                # (which are over the identity threshold): Hit gene's EC numbers and Best hit information
                li=raw_fp.readlines()
                best_hit_lines={}
                cpt=0

                for i, line in enumerate(li):
                    if re.match("^#",line):
                        pass
                    else:
                        li_split=line.split("\t")

                        if float(li_split[2]) >= id : # Parsing blast hit according to the identity threshold
                            cpt+=1
                            hit_gene_ids = li_split[1]
                            gene_uniprotKB_ID = hit_gene_ids.split('|')
                            gene_name = str(gene_uniprotKB_ID[1])
                            key = str(li_split[0])

                            if  key in best_hit_lines : # Give information about the best hit for each assessed gene -> Output dict.
                                hit_parsed = "{}\t{}\t{}".format(line.rstrip(),"SECONDARY_HITS",gene_ec[gene_name])
                            else:
                                best_hit_lines[key] = 1
                                hit_parsed = "{}\t{}\t{}".format(line.rstrip(),"BEST_HIT",gene_ec[gene_name])
                            
                            
                            filtered_fp.write(hit_parsed)
                        else:
                            pass
