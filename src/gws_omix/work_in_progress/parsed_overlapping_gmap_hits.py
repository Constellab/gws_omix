# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import task_decorator, File, FloatParam, IntParam, ConfigParams, TaskInputs, TaskOutputs, Utils
from ..base.omix_env_task import BaseOmixEnvTask
from ..file.gff3_file import GFF3File
from ..file.gtf_file import GTFFile

@task_decorator("ParsedOverlappingGmapHits", hide=True)
class ParsedOverlappingGmapHits(BaseOmixEnvTask):
    """
    ParsedOverlappingGmapHits class. Represents a process that wraps bedtools tools to select a gene feature when multiple are overlapping. To use with gmap output gff3 file.
    
    Configuration options
        * `A_overlap"`: Proportion of feature A overlapped by feature B (default=0.5, 0.0 would means no overlap). 
        * `B_overlap`: Proportion of feature B overlapped by feature A  (default=0.5, 0.0 would means no overlap).
    """
    input_specs = {
        'gmap_gff3_file': (GFF3File,)       
    }
    output_specs = {
        'Parsed_Gmap_gtf_file': (GTFFile,)
    }
    config_specs = {
        "A_overlap": FloatParam(default_value=0.5, min_value=0.0, max_value=1.0, short_description="Proportion of feature A overlapped by feature B (default=0.5, 0.0 would means no overlap)."),
        "B_overlap": FloatParam(default_value=0.5, min_value=0.0, max_value=1.0, short_description="Proportion of feature A overlapped by feature B (default=0.5, 0.0 would means no overlap)."),
    }   
    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        result_file = File()
        result_file.path = self._get_output_file_path(params)
        return {"Parsed_Gmap_gtf_file": result_file} 
   
    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        A_proportion = params["A_overlap"]
        B_proportion = params["B_overlap"]
        gmap_gff3 = params["gmap_gff3_file"]

        cmd = [
            "bedtools intersect  ",
            "-sorted -f ", A_proportion,
            " -F ", B_proportion,
            " -a <( sort -k1,1 -k4,4n ", gmap_gff3,
            "  | grep \"mRNA\") -b <( sort -k1,1 -k4,4n  ", gmap_gff3,
            "  | grep \"mRNA\") -wa -wb | awk '{if($9 != $18){ print $0}}   ", cov,
            #"perl -ne ' #keeps and return not best hits# compare $t[8] and  $t[19] for matches, cov and identity' | filterring out not best hit in the gff3 file  > ", tmp.parsed.gff3, 
            "gffread tmp.parsed.gff3  -T -o ", self._get_output_file_path(params),
            " ; rm tmp.parsed.gff3  ; "
        ]           
# add when overlapping feature on a given position, #1 : keep the gene feature with highest score then #2: if equals keep the largest gene feature
        return cmd

    def _get_output_file_path(self, params):
        return params["cdna_or_cds_fasta_file"]+ ".alligned_on." + params["uncompressed_genome_fasta_file"] + ".gmap_alignment.parsed.gtf"

#bedtools intersect -sorted -f 0.8 -F 0.8   -a <( sort -k1,1 -k4,4n ecoli_ealbertii.test.gff3  | grep "mRNA") -b <( sort -k1,1 -k4,4n ecoli_ealbertii.test.gff3 | grep "mRNA") -wa -wb | awk '{if($9 != $18){ print $0}}'