# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, IntParam, StrParam, FloatParam, ConfigParams, TaskInputs, TaskOutputs, Utils, Settings
from ..base.omix_env_task import BaseOmixEnvTask
from ..file.fasta_file import FastaFile
from ..file.blast_ec_file import BlastECFile

@task_decorator("BlastEC")
class BlastEC(BaseOmixEnvTask):
    """
    BlastEC class. Represents a process that wraps NCBI blast program. This version !!! ALLOWED !!! to get EC numbers for digital twins reconstruction.
    
    Configuration options
        * `taxo`: Kingdom name. Specify taxonomic groups to select a specific sub-set database (Faster) = bacteria, archaea, eukaryota, metazoa, chordata, mammalia, fungi, viridiplantae. [Default = all] = Slower.
        * `alignement_type`: Alignement type. Prot against Prot database (i.e blastp) or Translated Nucl against prot database (i.e blastx). Respectivly, options = PP, TNP. [Default: PP] ".
        * `num_alignments: Number of database sequences to show alignments for [Default: 10],
        * `evalue`: E-value to exclude results. Default = 0.00001 (i.e 1e-5).
        * `threads`: Multi threading options: number of threads to use (min=1, max=7). [Default =  4].
        * `idt`: Similarity/identity minimum percentage threshold to exclude results (min= 1, max= 100). [Default = 70].
        * `cov`: Coverage (see blast option -qcov_hsp_perc) minimum percentage threshold to exclude results (min= 1, max= 100). [Default = 70].

    """
   
    input_specs = {
        'fasta_file': (FastaFile,)
    }
    output_specs = {
        'filtered_blast_ec_file': (BlastECFile,)
    }
    config_specs = {
        "taxonomy": StrParam(allowed_values=["all", "prokaryota", "eukaryota", "animals", "fungi", "plant"],  short_description="Kingdom name : Specify kingdom to select the database (Faster) = prokaryota, eukaryota, animals, fungi or plant"),
        "alignment_type": StrParam(default_value="PP",allowed_values=["PP","TNP"], short_description="Type of alignement to perform : Prot against Prot database (i.e blastp) or Translated Nucl against prot database (i.e blastx). [Respectivly, options : PP, TNP ]. Default = PP"),
        "num_alignments": IntParam(default_value=10, min_value=1, max_value=250, short_description="Number of database sequences to show alignments for [Default: 10]"),
        "e_value": FloatParam(default_value=0.00001, min_value=0.0, short_description="E-value : Default = 0.00001 (i.e 1e-5)"),
        "threads": IntParam(default_value=4, min_value=2, short_description="Number of threads"),
        "idt": IntParam(default_value=70, min_value=1, max_value=100, short_description="Similarity/identity minimum percentage threshold to exclude results. [Default = 70]"),
        "cov": IntParam(default_value=70, min_value=1, max_value=100, short_description="Coverage (see blast option -qcov_hsp_perc) minimum percentage threshold to exclude results [Default = 70]"),
        "uniprot_db_dir": StrParam(default_value="", short_description="Location of the UniProtKB database") # TODO: set protected
    }
    
    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        # execute blast parsing command and ec number retrieving 
        taxo = params["taxonomy"]
        idt = params["idt"]
        local_uniprot_db_dir = params["uniprot_db_dir"]
        tab_file = os.path.join(local_uniprot_db_dir, taxo + ".uniprotKB.tab")        
        filtered_file_path = self._create_filtered_output_file(self._output_file_path, tab_file, idt)
        result_file = BlastECFile(path=filtered_file_path)
        return {"filtered_blast_ec_file": result_file}

    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        taxo  = params["taxonomy"]
        alignment  = params["alignment_type"]
        evalue  = params["e_value"]
        thread  = params["threads"]
        num_alignments  = params["num_alignments"]
        cov  = params["cov"]
        local_uniprot_db_dir = params["uniprot_db_dir"]

        fasta_file = inputs["fasta_file"]
        fasta_file_name = os.path.basename(fasta_file.path)

        datab_file_path = os.path.join(local_uniprot_db_dir, taxo + ".uniprotKB.faa")
        self._output_file_path = self._get_output_file_path(params["taxonomy"], fasta_file_name)

        script_file_dir = os.path.dirname(os.path.realpath(__file__))
        if alignment == "PP":
            cmd = [  
                "bash", 
                os.path.join(script_file_dir, "./sh/blastp_cmd.sh"),
                datab_file_path,
                fasta_file.path,
                evalue,
                thread,
                cov,
                num_alignments,
                self._output_file_path
            ]                 
        else:
            cmd = [  
                "bash", 
                os.path.join(script_file_dir, "./sh/blastx_cmd.sh"),
                datab_file_path,
                fasta_file.path,
                evalue,
                thread,
                cov,
                num_alignments,
                self._output_file_path
            ]  

        return cmd

    def _get_output_file_path(self, taxonomy, fasta_file_name) :
        return os.path.join(
            self.working_dir, 
            fasta_file_name + ".alligned_on." + taxonomy + ".blast_output"
        )

    def _create_filtered_output_file(self, blast_output_file, tabular_file, id)  :
        gene_ec={}
        hit_parsed={}

        with open(tabular_file, 'r') as lines: # Create dict. containing genes with their corresponding EC number(s)
            li=lines.readlines()
            for _, line in enumerate(li):
                if re.match("^#",line):
                    pass
                else:
                    li_split=line.split("\t")
                    if re.match("^$",str(li_split[7])):
                        gene_ec[li_split[0]]="NA\n"                   
                    else:
                        gene_ec[li_split[0]]=li_split[7]                    

#        for current_line in blast_output_file:

        filtered_file_path = blast_output_file + ".filtered.csv"

        with open(filtered_file_path, 'w+') as filtered_file_fp: 
            with open(blast_output_file, 'r') as raw_fp: 
                # Create dict. containing for each lines of the blast output 
                # (which are over the identity threshold): Hit gene's EC numbers and Best hit information
                li=raw_fp.readlines()
                best_hit_lines={}
                cpt=0

                for _, line in enumerate(li):
                    if re.match("^#",line):
                        pass
                    else:
                        li_split=line.split("\t")
                        if float(li_split[2]) >= float(id)  : # Parsing blast hit according to the identity threshold
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
                                                    
                            #filtered_dict[cpt]=hit_parsed
                            filtered_file_fp.write( hit_parsed )
                        else:
                            pass

        return filtered_file_path
