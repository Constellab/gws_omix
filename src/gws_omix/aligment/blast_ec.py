# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import re

from gws_core import (CondaShellProxy, ConfigParams, ConfigSpecs, FloatParam,
                      InputSpec, InputSpecs, IntParam, OutputSpec, OutputSpecs,
                      StrParam, Task, TaskFileDownloader, TaskInputs,
                      TaskOutputs, task_decorator)
from gws_omix.base_env.omix_env_task import BaseOmixEnvHelper

from ..file.blast_ec_file import BlastECFile
from ..file.fasta_file import FastaFile

# from ..utils._settings import Settings


@task_decorator("BlastEC", human_name="Blast to EC-number Annotator",
                short_description="BlastEC is an homology based EC-number (Enzyme id) annotator which used Blastp/x and UniProtKB-db.")
class BlastEC(Task):
    """
    BlastEC class.

    Represents a process that wraps NCBI blast program. This version !!! ALLOWED !!! to get EC numbers for digital twins reconstruction.

    Configuration options
        * `taxonomy`: Specify the tax group to select the dedicated database.
        * `alignement_type`: Alignement type. Prot against Prot database (i.e blastp) or Translated Nucl against prot database (i.e blastx). Respectivly, options = PP, TNP. [Default: PP] ".
        * `num_alignments: Number of database sequences to show alignments for [Default: 10],
        * `evalue`: E-value to exclude results. Default = 0.00001 (i.e 1e-5).
        * `threads`: Multi threading options: number of threads to use (min=1, max=7). [Default =  1].
        * `idt`: Similarity/identity minimum percentage threshold to exclude results (min= 1, max= 100). [Default = 70].
        * `cov`: Coverage (see blast option -qcov_hsp_perc) minimum percentage threshold to exclude results (min= 1, max= 100). [Default = 70].

    """

    DB_LOCATION = 'https://storage.gra.cloud.ovh.net/v1/AUTH_a0286631d7b24afba3f3cdebed2992aa/opendata/omix/db'
    TAX_DICT = {
        "archaea": "uniprot-taxonomy_v1_2157",
        "bacteria": "uniprot-taxonomy_v1_2",
        "chordata": "uniprot-taxonomy_v1_7711",
        "eukaryota": "uniprot-taxonomy_v1_2759",
        "fungi": "uniprot-taxonomy_v1_4751",
        "mammalia": "uniprot-taxonomy_v1_40674",
        "metazoa": "uniprot-taxonomy_v1_33208",
        "viridiplantae": "uniprot-taxonomy_v1_33090",
        "virus": "uniprot-taxonomy_v1_10239"
    }

    OUT_FMT = '7 qaccver saccver pident qcovs qcovhsp length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore'

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
        "threads": IntParam(default_value=1, min_value=1, short_description="Number of threads"),
        "idt": IntParam(default_value=70, min_value=1, max_value=100, short_description="Similarity/identity minimum percentage threshold to exclude results. [Default = 70]"),
        "cov": IntParam(default_value=70, min_value=1, max_value=100, short_description="Coverage (see blast option -qcov_hsp_perc) minimum percentage threshold to exclude results [Default = 70]"),
        # TODO: set protected
        # "uniprot_db_dir": StrParam(default_value="", short_description="Location of the UniProtKB database")
    }

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        shell_proxy: CondaShellProxy = BaseOmixEnvHelper.create_proxy(
            self.message_dispatcher)

        taxo = params["taxonomy"]
        alignment = params["alignment_type"]
        evalue = params["e_value"]
        thread = params["threads"]
        num_alignments = params["num_alignments"]
        cov = params["cov"]
        # local_uniprot_db_dir = params["uniprot_db_dir"]

        # DB location
        file_prefix = self.TAX_DICT[taxo]

        fasta_file_path = self._download_fasta(file_prefix)
        blastdb_folder_path = self._download_blastdb(file_prefix)

        # add fasta file path, and blastdb_folder_path to working dir using symlink
        # because blast command need to have the fasta file and the blastdb in the same folder
        shell_proxy.run(["ln", "-s", fasta_file_path, './'])
        shell_proxy.run(["ln", "-s", blastdb_folder_path + '/*', './'])

        fasta_file: FastaFile = inputs["fasta_file"]
        fasta_file_name = os.path.basename(fasta_file.path)

        output_file_path = os.path.join(
            shell_proxy.working_dir,
            fasta_file_name + ".alligned_on." + taxo + ".blast_output")

        if alignment == "PP":
            cmd = [
                "blastp",
                "-db", "uniprot*.fasta",
                "-query", fasta_file_name,
                "-evalue", evalue,
                "-num_threads", thread,
                "-qcov_hsp_perc", cov,
                "-num_alignments", num_alignments,
                "-outfmt", self.OUT_FMT,
                "show_gis",
                "-task", "balstp-fast",
                "-out", output_file_path
            ]
        else:
            cmd = [
                "blastx",
                "-db", "uniprot*.fasta",
                "-query", fasta_file_name,
                "-evalue", evalue,
                "-num_threads", thread,
                "-qcov_hsp_perc", cov,
                "-num_alignments", num_alignments,
                "-outfmt", self.OUT_FMT,
                "show_gis",
                "-task", "balstp-fast",
                "-out", output_file_path
            ]

        # call the command
        shell_proxy.run(cmd, shell_mode=True)

        # execute blast parsing command and ec number retrieving
        idt = params["idt"]

        tab_file_path = self._download_tab(file_prefix)
        filtered_file_path = self._create_filtered_output_file(
            output_file_path, tab_file_path, idt)
        result_file = BlastECFile(path=filtered_file_path)

        return {"filtered_blast_ec_file": result_file}

    def _download_fasta(self, file_prefix: str) -> str:
        file_downloader = TaskFileDownloader(
            BlastEC.get_brick_name(), self.message_dispatcher)

        # Gest fasta
        self.log_info_message("Downloading fasta file")
        fasta_file_name = file_prefix + '.fasta.gz'
        db_location = os.path.join(self.DB_LOCATION)
        return file_downloader.download_file_if_missing(
            db_location, fasta_file_name, decompress_file=True)

    def _download_tab(self, file_prefix: str) -> str:
        file_downloader = TaskFileDownloader(
            BlastEC.get_brick_name(), self.message_dispatcher)

        # Get tab
        self.log_info_message("Downloading tab file")
        tab_file_name = file_prefix + '.tab.gz'
        db_location = os.path.join(self.DB_LOCATION, tab_file_name)
        return file_downloader.download_file_if_missing(
            db_location, tab_file_name, decompress_file=True)

    def _download_blastdb(self, file_prefix: str) -> str:
        file_downloader = TaskFileDownloader(
            BlastEC.get_brick_name(), self.message_dispatcher)

        # Get blastdb
        self.log_info_message("Downloading blastdb file")
        blastdb_file_name = file_prefix + '.fasta_blast_index.tar.gz'
        db_location = os.path.join(self.DB_LOCATION, blastdb_file_name)
        return file_downloader.download_file_if_missing(
            db_location, blastdb_file_name, decompress_file=True)

    def _create_filtered_output_file(self, blast_output_file: str, tabular_file: str, id):
        gene_ec = {}
        hit_parsed = {}

        # Create dict. containing genes with their corresponding EC number(s)
        with open(tabular_file, 'r') as lines:
            for line in lines:
                if re.match("^#", line):
                    pass
                else:
                    li_split = line.split("\t")
                    if re.match("^$", str(li_split[7])):
                        gene_ec[li_split[0]] = "NA\n"
                    else:
                        gene_ec[li_split[0]] = li_split[7] + \
                            "\n"  # ! : missing+ "\n"

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
                        # Parsing blast hit according to the identity threshold
                        if float(li_split[2]) >= float(id):
                            cpt += 1
                            hit_gene_ids = li_split[1]
                            # gene_uniprotKB_ID = hit_gene_ids  # hit_gene_ids.split('|')
                            # str(gene_uniprotKB_ID[1])
                            gene_name = str(hit_gene_ids)
                            key = str(li_split[0])
                            # Give information about the best hit for each assessed gene -> Output dict.
                            if key in best_hit_lines:
                                hit_parsed = "{}\t{}\t{}".format(
                                    line.rstrip(), "SECONDARY_HITS", gene_ec[gene_name])
                            else:
                                best_hit_lines[key] = 1
                                hit_parsed = "{}\t{}\t{}".format(
                                    line.rstrip(), "BEST_HIT", gene_ec[gene_name])

                            # filtered_dict[cpt]=hit_parsed
                            filtered_file_fp.write(hit_parsed)
                        else:
                            pass

        return filtered_file_path
