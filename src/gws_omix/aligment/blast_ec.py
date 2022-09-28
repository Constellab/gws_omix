# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import csv
import json
import os
import re

from gws_core import (File, FloatParam, InputSpec, IntParam, OutputSpec,
                      Settings, StrParam, TaskInputs, TaskOutputs, Utils,
                      task_decorator)
from gws_core.config.config_types import ConfigParams, ConfigSpecs
from gws_core.io.io_spec import InputSpec, OutputSpec
from gws_core.io.io_spec_helper import InputSpecs, OutputSpecs

from ..base_env.omix_env_task import BaseOmixEnvTask
from ..file.blast_ec_file import BlastECFile
from ..file.fasta_file import FastaFile


@task_decorator("BlastEC")
class BlastEC(BaseOmixEnvTask):
    """
    BlastEC class.

    Represents a process that wraps NCBI blast program. This version !!! ALLOWED !!! to get EC numbers for digital twins reconstruction.

    Configuration options
        * `taxo`: Kingdom name. Specify taxonomic groups to select a specific sub-set database (Faster) = bacteria, archaea, eukaryota, metazoa, chordata, mammalia, fungi, viridiplantae. [Default = all] = Slower.
        * `alignement_type`: Alignement type. Prot against Prot database (i.e blastp) or Translated Nucl against prot database (i.e blastx). Respectivly, options = PP, TNP. [Default: PP] ".
        * `num_alignments: Number of database sequences to show alignments for [Default: 10],
        * `evalue`: E-value to exclude results. Default = 0.00001 (i.e 1e-5).
        * `threads`: Multi threading options: number of threads to use (min=1, max=7). [Default =  4].
        * `idt`: Similarity/identity minimum percentage threshold to exclude results (min= 1, max= 100). [Default = 70].
        * `cov`: Coverage (see blast option -qcov_hsp_perc) minimum percentage threshold to exclude results (min= 1, max= 100). [Default = 70].

    """
    TAX_DICT = {
        "archaea": "/data/gws_omix/opendata/uniprot-taxonomy_v1_2157",
        "bacteria": "/data/gws_omix/opendata/uniprot-taxonomy_v1_2",
        "chordata": "/data/gws_omix/opendata/uniprot-taxonomy_v1_7711",
        "eukaryota": "/data/gws_omix/opendata/uniprot-taxonomy_v1_2759",
        "fungi": "/data/gws_omix/opendata/uniprot-taxonomy_v1_4751",
        "mammalia": "/data/gws_omix/opendata/uniprot-taxonomy_v1_40674",
        "metazoa": "/data/gws_omix/opendata/uniprot-taxonomy_v1_33208",
        "viridiplantae": "/data/gws_omix/opendata/uniprot-taxonomy_v1_33090",
        "virus": "/data/gws_omix/opendata/uniprot-taxonomy_v1_10239"
    }
    input_specs: InputSpecs = {
        'fasta_file': InputSpec(FastaFile, human_name="Fasta File", short_description="Fasta Input File")
    }
    output_specs: OutputSpecs = {
        'filtered_blast_ec_file': OutputSpec(BlastECFile, human_name="Blast results", short_description="Blast results")
    }
    config_specs: ConfigSpecs = {
        "taxonomy": StrParam(allowed_values=["fungi"],  short_description="Specify the tax group to select the dedicated database"),
        # "all", "prokaryota", "eukaryota", "animals", "fungi", "plant"
        "alignment_type": StrParam(default_value="PP", allowed_values=["PP", "TNP"], short_description="Type of alignement to perform : Prot against Prot database (i.e blastp) or Translated Nucl against prot database (i.e blastx). [Respectivly, options : PP, TNP ]. Default = PP"),
        "num_alignments": IntParam(default_value=10, min_value=1, max_value=250, short_description="Number of database sequences to show alignments for [Default: 10]"),
        "e_value": FloatParam(default_value=0.00001, min_value=0.0, short_description="E-value : Default = 0.00001 (i.e 1e-5)"),
        "threads": IntParam(default_value=4, min_value=2, short_description="Number of threads"),
        "idt": IntParam(default_value=70, min_value=1, max_value=100, short_description="Similarity/identity minimum percentage threshold to exclude results. [Default = 70]"),
        "cov": IntParam(default_value=70, min_value=1, max_value=100, short_description="Coverage (see blast option -qcov_hsp_perc) minimum percentage threshold to exclude results [Default = 70]"),
        # TODO: set protected
        # "uniprot_db_dir": StrParam(default_value="", short_description="Location of the UniProtKB database")
    }

    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        # execute blast parsing command and ec number retrieving
        taxo = params["taxonomy"]
        idt = params["idt"]
        #local_uniprot_db_dir = params["uniprot_db_dir"]

        tab_file = os.path.join(self.TAX_DICT[taxo] + ".tab")
        filtered_file_path = self._create_filtered_output_file(self.output_file_path, tab_file, idt)
        result_file = BlastECFile(path=filtered_file_path)
        return {"filtered_blast_ec_file": result_file}

    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        taxo = params["taxonomy"]
        alignment = params["alignment_type"]
        evalue = params["e_value"]
        thread = params["threads"]
        num_alignments = params["num_alignments"]
        cov = params["cov"]
        #local_uniprot_db_dir = params["uniprot_db_dir"]
        fasta_file = inputs["fasta_file"]
        fasta_file_name = os.path.basename(fasta_file.path)

        #datab_file_path = os.path.join(local_uniprot_db_dir, taxo + ".uniprotKB.faa")
        # datab_file_path = os.path.join(
        #     self.TAX_DICT[taxo] + ".fasta_blast_index", self.TAX_DICT[taxo] + ".fasta")
        datab_file_path = os.path.join(self.TAX_DICT[taxo] + ".fasta")
        self.output_file_path = self._get_output_file_path(taxo, fasta_file_name)

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
                self.output_file_path
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
                self.output_file_path
            ]

        return cmd

    def _get_output_file_path(self, taxonomy, fasta_file_name):
        return os.path.join(
            self.working_dir,
            fasta_file_name + ".alligned_on." + taxonomy + ".blast_output"
        )

    def _create_filtered_output_file(self, blast_output_file, tabular_file, id):
        gene_ec = {}
        hit_parsed = {}

        with open(tabular_file, 'r') as lines:  # Create dict. containing genes with their corresponding EC number(s)
            li = lines.readlines()
            for _, line in enumerate(li):
                if re.match("^#", line):
                    pass
                else:
                    li_split = line.split("\t")
                    if re.match("^$", str(li_split[7])):
                        gene_ec[li_split[0]] = "NA\n"
                    else:
                        gene_ec[li_split[0]] = li_split[7]

        filtered_file_path = blast_output_file + ".filtered.csv"

        with open(filtered_file_path, 'w+') as filtered_file_fp:
            with open(blast_output_file, 'r') as raw_fp:
                # Create dict. containing for each lines of the blast output
                # (which are over the identity threshold): Hit gene's EC numbers and Best hit information
                li = raw_fp.readlines()
                best_hit_lines = {}
                cpt = 0

                for _, line in enumerate(li):
                    if re.match("^#", line):
                        pass
                    else:
                        li_split = line.split("\t")
                        if float(li_split[2]) >= float(id):  # Parsing blast hit according to the identity threshold
                            cpt += 1
                            hit_gene_ids = li_split[1]
                            gene_uniprotKB_ID = hit_gene_ids.split('|')
                            gene_name = str(gene_uniprotKB_ID[1])
                            key = str(li_split[0])
                            # Give information about the best hit for each assessed gene -> Output dict.
                            if key in best_hit_lines:
                                hit_parsed = "{}\t{}\t{}".format(line.rstrip(), "SECONDARY_HITS", gene_ec[gene_name])
                            else:
                                best_hit_lines[key] = 1
                                hit_parsed = "{}\t{}\t{}".format(line.rstrip(), "BEST_HIT", gene_ec[gene_name])

                            # filtered_dict[cpt]=hit_parsed
                            filtered_file_fp.write(hit_parsed)
                        else:
                            pass

        return filtered_file_path
