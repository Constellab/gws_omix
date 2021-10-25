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

@task_decorator("Blast")
class Blast(BaseOmixEnvTask):
    """
    Blast class. Represents a process that wraps NCBI blast program. This version !!! DOES NOT ALLOWED !!! to get EC numbers for digital twins reconstruction.
    
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
        'filtered_blast_file': (File,)
    }
    config_specs = {
        "taxo": StrParam(allowed_values=["all", "prokaryota", "eukaryota", "animals", "fungi", "plant"],  short_description="Kingdom name : Specify kingdom to select the database (Faster) = prokaryota, eukaryota, animals, fungi or plant"),
        "alignment_type": StrParam(default_value="PP",allowed_values=["PP","TNP"], short_description="Type of alignement to perform : Prot against Prot database (i.e blastp) or Translated Nucl against prot database (i.e blastx). [Respectivly, options : PP, TNP ]. Default = PP"),
        "num_alignments": IntParam(default_value=10, min_value=1, max_value=250, short_description="Number of database sequences to show alignments for [Default: 10]"),
        "e_value": FloatParam(default_value=0.00001, min_value=0.0, short_description="E-value : Default = 0.00001 (i.e 1e-5)"),
        "threads": IntParam(default_value=4, min_value=2, short_description="Number of threads"),
        "idt": IntParam(default_value=70, min_value=1, max_value=100, short_description="Similarity/identity minimum percentage threshold to exclude results. [Default = 70]"),
        "cov": IntParam(default_value=70, min_value=1, max_value=100, short_description="Coverage (see blast option -qcov_hsp_perc) minimum percentage threshold to exclude results [Default = 70]")
    }

    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        # execute blast cmd
        result_file = File()
        result_file.path = self._get_output_file_path(params)
        
        # execute blast parsing command and ec number retrieving 
        idt  = params["idt"]       
        
        result_file_2 = File()
        result_file_2.path = self._create_filtered_output_file(result_file.path, idt)

        return {"filtered_blast_file": result_file_2}
  
    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        taxo  = params["taxo"]
        alignment  = params["alignment_type"]
        evalue  = params["e_value"]
        thread  = params["threads"]
        num_alignments  = params["num_alignments"]
        cov  = params["cov"]
        fasta_file = params["fasta_file"]
        datab_dir = "$blast_index_dir"
        datab_file_path = os.path.join(datab_dir, taxo + ".uniprotKB.faa")

        if alignment == "PP":
            blast_type = "blastp"
        else:
            blast_type = "blastx"

        cmd = [ 
            blast_type, 
            " -db ", datab_file_path,        
            " -query ", fasta_file.path,
            " -evalue ", evalue,
            " -num_threads ", thread,
            " -qcov_hsp_perc ", cov,
            " -num_alignments ", num_alignments,
            "-outfmt  \"7 qaccver saccver pident qcovs qcovhsp length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore\" -show_gis -task \"blastp-fast\" -out ", self._get_output_file_path(params)
        ]
        return cmd

    def _get_output_file_path(self, params):
        return params["fasta_file"]+ ".alligned_on." + params["taxo"] + ".blast_output.csv"

    def _create_filtered_output_file(self, blast_output_file, id):
        hit_parsed={}
        filtered_fp={}
        with open(blast_output_file, 'r') as raw_fp: 
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
                        key = str(li_split[0])

                        if  key in best_hit_lines : # Give information about the best hit for each assessed gene -> Output dict.
                            hit_parsed = "{}\t{}".format(line.rstrip(),"SECONDARY_HITS")
                        else:
                            best_hit_lines[key] = 1
                            hit_parsed = "{}\t{}".format(line.rstrip(),"BEST_HIT")

                        filtered_fp[cpt]=hit_parsed
                    else:
                        pass

        return filtered_fp
