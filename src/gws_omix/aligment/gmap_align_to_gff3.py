# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import json
import os
import re
import csv

from gws_core import (task_decorator, File, IntParam, StrParam, FloatParam, 
                        ConfigParams, TaskInputs, TaskOutputs, Utils, Folder, OptionalIn, BoolParam)
from ..base_env.omix_env_task import BaseOmixEnvTask
from ..file.fasta_file import FastaFile
from ..file.gff3_file import GFF3File


@task_decorator("GmapAlignGFF3")
class GmapAlignGFF3(BaseOmixEnvTask):
    """
    GmapAlignGFF3 class. 
    
    Represents a task that wraps Gmap aligment tool. Gmap index is Mandatory to use Gmap tools.
    
    Configuration options:
        * `threads`: Multi threading options: number of threads to use. [Default =  8]
        * `min-identity"`: Do not print alignments with identity less this value (default=70.0, 0.0 would means no filtering). Note that chimeric alignments will be output regardless of this filter. 
        * `min-trimmed-coverage`: Do not print alignments with trimmed coverage less than this value (default=70.0, 0.0 would means no filtering). Note that chimeric alignments will be output regardless of this filter. 
        * `max-hit-number`: Maximum number of hits to show (default 5).  If set to 1, GMAP will not report chimeric alignments, since those imply two hits. If you want a single alignment plus chimeric alignments, then set this to be 0. [Default =  5] 
        * `cross-species`: Use a more sensitive search for canonical splicing, which helps especially for cross-species alignments and other difficult cases (genome for a far-related species/family...). [Default =  No]
    """
    input_specs = {
        'uncompressed_genome_fasta_file': (FastaFile,),
#        'uncompressed_genome_fasta_file_index_dir': OptionalIn(Folder,),
        'cdna_or_cds_fasta_file': (FastaFile,)          
    }
    output_specs = {
        'gmap_gff3_file': (GFF3File,)
    }
    config_specs = {
        "threads": IntParam(default_value=8, min_value=1, short_description="Number of threads [Default =  8] "),
        "min-identity": FloatParam(default_value=0.7, max_value=1.0, min_value=0.0, short_description="Do not print alignments with identity less this value (default=0.7, 0.0 would means no filtering). Note that chimeric alignments will be output regardless of this filter. "),
        "min-trimmed-coverage": FloatParam(default_value=0.7, min_value=0.0, max_value=1.0, short_description="Do not print alignments with trimmed coverage less than this value (default=0.7, 0.0 would means no filtering). Note that chimeric alignments will be output regardless of this filter. "),
        "max-hit-number": IntParam(default_value=5, min_value=0, short_description="Maximum number of hits to show (default 5).  If set to 1, GMAP will not report chimeric alignments, since those imply two hits. If you want a single alignment plus chimeric alignments, then set this to be 0. [Default =  5] "),
        "cross-species": StrParam(default_value="No", allowed_values=["Yes","No"], short_description="Use a more sensitive search for canonical splicing, which helps especially for cross-species alignments and other difficult cases (genome for a far-related species/family...). [Default =  No]"),
        "alt-start-codons": StrParam(default_value="No", allowed_values=["Yes","No"], short_description="Also, use the alternate initiation codons (see http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). By default, without this option, only ATG is considered an initiation codon [Default =  No]"),
        "fulllength": StrParam(default_value="No", allowed_values=["Yes","No"], short_description="Assume full-length protein, starting with Met (ATG codon). [Default =  No]"),
        "use_default_index_dirs": BoolParam(default_value=True, visibility="protected", short_description="True to use the default index dirs provided by us. False to use the index provided as input. In any case, if an index is provided in inputs, it is used.")
    }
   
    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        result_file = GFF3File()
        result_file.path = self._output_file_path
        return {"gmap_gff3_file": result_file} 
   
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
        genome_fasta_dir = inputs["uncompressed_genome_fasta_file"]

        fa_file = inputs["cdna_or_cds_fasta_file"]
        fasta_file_name = os.path.basename(fa_file.path)
        self._output_file_path = self._get_output_file_path(fasta_file_name, genome_fasta_file_name)
        script_file_dir = os.path.dirname(os.path.realpath(__file__))

        cmd = [
            "bash", 
            os.path.join(script_file_dir, "./sh/gmap_gff3_cmd.sh"),
            thread,
            idt,
            cov,
            hit_nbr,
            crs_species,
            alt_start,
            full_lgth,
            genome_fasta.path,
            fa_file.path,
            self._output_file_path,
            genome_fasta_dir.path
        ]
        return cmd

    def _get_output_file_path(self, fasta, genome):
        return os.path.join(
            self.working_dir, 
            fasta + ".alligned_on." + genome + ".gmap_alignment.gff3"
        )
