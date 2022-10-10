# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import csv
import json
import os
import re

from gws_core import (File, FloatParam, InputSpec, IntParam, OutputSpec,
                      Settings, StrParam, Table, TaskInputs, TaskOutputs,
                      Utils, task_decorator)
from gws_core.config.config_types import ConfigParams, ConfigSpecs
from gws_core.io.io_spec import InputSpec, OutputSpec
from gws_core.io.io_spec_helper import InputSpecs, OutputSpecs

from ..base_env.interproscan_env_task import InterProScanEnvTask
from ..file.blast_ec_file import BlastECFile
from ..file.fasta_file import FastaFile
from ..utils._requests import Requests

# from ..utils._settings import Settings


@task_decorator("InterProScan", human_name="InterProScan Annotator",
                short_description="InterProScan is an homology based protein annotator (see https://www.ebi.ac.uk/interpro/search/sequence/ )")
class InterProScan(InterProScanEnvTask):
    """
    InterProScan class.

    """
    input_specs: InputSpecs = {
        'fasta_file': InputSpec(FastaFile, human_name="Fasta File", short_description="Fasta Input File")
    }
    output_specs: OutputSpecs = {
        'output_table': OutputSpec(Table, human_name="Interproscan_results", short_description="Interproscan results")
    }
    config_specs: ConfigSpecs = {
        "taxonomy": StrParam(allowed_values=["fungi"],  short_description="Specify the tax group to select the dedicated database"),
        # "all", "prokaryota", "eukaryota", "animals", "fungi", "plant"
        "alignment_type": StrParam(default_value="PP", allowed_values=["PP", "TNP"], short_description="Type of alignement to perform : Prot against Prot database (i.e blastp) or Translated Nucl against prot database (i.e blastx). [Respectivly, options : PP, TNP ]. Default = PP"),
        "num_alignments": IntParam(default_value=10, min_value=1, max_value=250, short_description="Number of database sequences to show alignments for [Default: 10]"),
        "e_value": FloatParam(default_value=0.00001, min_value=0.0, short_description="E-value : Default = 0.00001 (i.e 1e-5)"),
        "threads": IntParam(default_value=1, min_value=1, short_description="Number of threads"),
        "idt": IntParam(default_value=70, min_value=1, max_value=100, short_description="Similarity/identity minimum percentage threshold to exclude results. [Default = 70]"),
        "cov": IntParam(default_value=70, min_value=1, max_value=100, short_description="Coverage (see blast option -qcov_hsp_perc) minimum percentage threshold to exclude results [Default = 70]"),
        # TODO: set protected
        # "uniprot_db_dir": StrParam(default_value="", short_description="Location of the UniProtKB database")
    }

    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        # execute blast parsing command and ec number retrieving
        taxo = params["taxonomy"]
        idt = params["idt"]
        fasta_file = inputs["fasta_file"]
        fasta_file_name = os.path.basename(fasta_file.path)
        # local_uniprot_db_dir = params["uniprot_db_dir"]

        tab_file_path = os.path.join(self.TAX_DICT[taxo] + ".tab")
        output_file_path = os.path.join(
            self.working_dir,
            fasta_file_name + ".alligned_on." + taxo + ".blast_output")
        filtered_file_path = self._create_filtered_output_file(output_file_path, tab_file_path, idt)
        result_file = BlastECFile(path=filtered_file_path)
        return {"filtered_blast_ec_file": result_file}

    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        settings = Settings.retrieve()
        taxo = params["taxonomy"]
        alignment = params["alignment_type"]
        evalue = params["e_value"]
        thread = params["threads"]
        num_alignments = params["num_alignments"]
        cov = params["cov"]
        # local_uniprot_db_dir = params["uniprot_db_dir"]

        # Check for DB existance
        test_filename = "gws_omix:uniprot_" + taxo + "_fasta_file"
        test_filename_2 = "gws_omix:uniprot_" + taxo + "_tab_file"
        test_filename_3 = "gws_omix:uniprot_" + taxo + "_blastdb_file"
        # test_path = settings["variables"][test_filename]
        test_path = settings.get_variable(test_filename)
        test_path = test_path.replace('.fasta.gz', '.fasta')
        test_path_2 = settings.get_variable(test_filename_2)
        test_path_2 = test_path_2.replace('.tab.gz', '.tab')
        test_path_3 = settings.get_variable(test_filename_3)
        test_path_3 = test_path_3.replace('.fasta_blast_index.tar.gz', '.fasta_blast_index')

        if os.path.exists(test_path) and os.path.exists(test_path_2) and os.path.exists(test_path_3):  # Pull BLASTdb
            print("# This DB already exists : start performing BLASTEC")
        else:
            print("# This DB does not exists : start downloading")
            # UniProtKB taxa_level selected DB

            # Get fasta
            taxa_prefix = "gws_omix:uniprot_" + taxo
            url_taxa_name = taxa_prefix + "_fasta_url"
            file_taxa_name = taxa_prefix + "_fasta_file"
            url = settings.get_variable(url_taxa_name)
            dest_path = settings.get_variable(file_taxa_name)
            Requests.download(url, dest_path)

            # Get blast DB
            taxa_prefix = "gws_omix:uniprot_" + taxo
            url_taxa_name = taxa_prefix + "_blastdb_url"
            file_taxa_name = taxa_prefix + "_blastdb_file"
            url = settings.get_variable(url_taxa_name)
            dest_path = settings.get_variable(file_taxa_name)
            Requests.download(url, dest_path)

            # Get EC tab file
            taxa_prefix = "gws_omix:uniprot_" + taxo
            url_taxa_name = taxa_prefix + "_tab_url"
            file_taxa_name = taxa_prefix + "_tab_file"
            url = settings.get_variable(url_taxa_name)
            dest_path = settings.get_variable(file_taxa_name)
            Requests.download(url, dest_path)

        fasta_file = inputs["fasta_file"]
        fasta_file_name = os.path.basename(fasta_file.path)

        datab_prefix_path = os.path.join(self.TAX_DICT[taxo])
        output_file_path = os.path.join(
            self.working_dir,
            fasta_file_name + ".alligned_on." + taxo + ".blast_output")

        script_file_dir = os.path.dirname(os.path.realpath(__file__))
        if alignment == "PP":
            cmd = [
                "bash",
                os.path.join(script_file_dir, "./sh/blastp_cmd.sh"),
                datab_prefix_path,
                fasta_file.path,
                evalue,
                thread,
                cov,
                num_alignments,
                output_file_path,
                self.OPENDATA_DIR
            ]
        else:
            cmd = [
                "bash",
                os.path.join(script_file_dir, "./sh/blastx_cmd.sh"),
                datab_prefix_path,
                fasta_file.path,
                evalue,
                thread,
                cov,
                num_alignments,
                output_file_path,
                self.OPENDATA_DIR
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
            for line in lines:
                if re.match("^#", line):
                    pass
                else:
                    li_split = line.split("\t")
                    if re.match("^$", str(li_split[7])):
                        gene_ec[li_split[0]] = "NA\n"
                    else:
                        gene_ec[li_split[0]] = li_split[7] + "\n"  # ! : missing+ "\n"

        filtered_file_path = blast_output_file + ".filtered.csv"

        with open(filtered_file_path, 'w+') as filtered_file_fp:
            with open(blast_output_file, 'r') as raw_fp:
                # Create dict. containing for each lines of the blast output
                # (which are over the identity threshold): Hit gene's EC numbers and Best hit information
                # li = raw_fp.readlines()
                best_hit_lines = {}
                cpt = 0
                for line in raw_fp:
                    # for _, line in enumerate(li):
                    if re.match("^#", line):
                        pass
                    else:
                        li_split = line.split("\t")
                        if float(li_split[2]) >= float(id):  # Parsing blast hit according to the identity threshold
                            cpt += 1
                            hit_gene_ids = li_split[1]
                            # gene_uniprotKB_ID = hit_gene_ids  # hit_gene_ids.split('|')
                            gene_name = str(hit_gene_ids)  # str(gene_uniprotKB_ID[1])
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
