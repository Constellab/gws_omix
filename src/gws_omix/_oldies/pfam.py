# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core.cond import JSONData
from gws.typing import Path
from gws.shell import Shell
from gws.file import File
from gws.logger import Error

from .data import PfamFile



#TO DO: pfam without ec number retrieving
#class Pfam(Shell):

class Pfam_ec(Shell):
    """
    Pfam class. Represents a process that wraps PFAM program.
    
    Configuration options
        * `threads`: Multi threading options: number of threads to use (min=1, max=7). [Default =  4].

    """

    input_specs = { 'fasta_file': (File,) }
    output_specs = { 'pfam_file': (File,) }
    config_specs = {
        "fasta": {"type": str, "Description": "Fasta input file" },
        "threads": {"type": int, "default": 4, "min": 1, "max": 7, "Description": "Number of threads"}
    }

    PFAM_DOMAIN_COL_INDEX = 0
    EC_NUMBER_COL_INDEX = 1   

    
    def after_command(self, stdout: str=None): 
        self.output["pfam_file"] = File( path=self.output_file ) 
        self.output_file_EC  = get_EC_number(self, self.output_file)
        self.output["pfam_file_EC"] = File( path=self.output_file_EC )        
        self.data["stdout"] = stdout
    
    def build_command(self) -> list:
        data_dir = Settings.get_data_dir()        
        brick_dir = self.get_brick_dir("omix")    
        bin_file = os.path.join(brick_dir, "bin", "pfam" , "pfam_scan.pl")
        datab_dir = os.path.join(data_dir ,"prod","omix","pfam")    
        fasta = self.get_param("fasta")
        thread = self.get_param("threads")


        self.output_file = os.path.join(self.cwd.name,".pfam.csv")
        
        cmd = [ 
        bin_file,
        "-dir ", datab_dir,
        "-fasta", fasta, 
        "-cpu ", thread,
        " > ", self.output_file
        ]
        
        return cmd


    def get_EC_number(self, pfam_output):
        data_dir = Settings.get_data_dir() 
        datab_domain_to_ec_number = os.path.join(data_dir ,"prod","omix","pfam","EC-Pfam_calculated_associations_Extended.csv")
        pfam_ec={}
        annotation_parsed={}

        with open(datab_domain_to_ec_number, 'r') as lines: # Create dict. containing pfam domain with their corresponding EC number
            li=lines.readlines()
            for line in li:                
                li_split=line.split("\t")
                pfam_ec[li_split[PFAM_DOMAIN_COL_INDEX]]=li_split[EC_NUMBER_COL_INDEX]                    

        with open(pfam_output, 'r') as lines: # Create dict. containing for each lines of the pfam annotation output file : gene <-> pfam domains(s) <-> EC number(s)
            li=lines.readlines()
            
            regexp_pfam = re.compile(r'.*(PF\d+)\..*')
            regexp_gene = re.compile(r'^([^\s]+)\s+.*$')
            for i, line in enumerate(li):  
                if re.match("^#",line):
                    pass
                elif re.match("^$",line):
                    pass
                else:                          
                    pfam_domain_id = regexp_pfam.findall(line)
                    gene_name = regexp_gene.findall(line)
                    gene = gene_name[0]
                    domain = pfam_domain_id[0]              
                    ec_nb = pfam_ec[domain]

                    if domain in pfam_ec:
                        annotation_parsed[i]= "{}\t{}\t{}".format(gene,domain,ec_nb)
                    else:
                        annotation_parsed[i]= "{}\t{}\t{}".format(gene,domain,"NA")                                

        return annotation_parsed

