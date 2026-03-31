# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

import pandas as pd

from gws_core import (
    ConfigParams,
    ConfigSpecs,
    File,
    FloatParam,
    InputSpec,
    InputSpecs,
    IntParam,
    OutputSpec,
    OutputSpecs,
    ShellProxy,
    StrParam,
    Table,
    TableImporter,
    Task,
    TaskFileDownloader,
    TaskInputs,
    TaskOutputs,
    task_decorator,
)


@task_decorator("Diamond", human_name="Diamond",
                short_description="Accelerated BLAST compatible local sequence aligner.")
class Diamond(Task):
    """
    DIAMOND, an evolutionary leap in bioinformatics, presenting an enhanced version that transcends the boundaries of its predecessors while matching the sensitivity of the gold standard BLASTP.
    It uses UniProtKB/Swiss-Prot database.

    Configuration options
        * `input_type_value`: Alignement type. Prot against Prot database (i.e blastp) or Translated Nucl against prot database (i.e blastx). Respectivly, options = prot, nuc. [Default: nuc] ".
        * `evalue_value`: E-value to exclude results. Default = 0.00001 (i.e 1e-5).
        * `num_threads`: Multi threading options: number of threads to use (min=1, max=7). [Default =  1].
        * `query_cover_value: Coverage (see blast option -qcov_hsp_perc) minimum percentage threshold to exclude results (min= 1, max= 100). [Default = 70]

    """

    DB_LOCATION = 'https://storage.gra.cloud.ovh.net/v1/AUTH_a0286631d7b24afba3f3cdebed2992aa/opendata/omix/db/swissprot_database.zip'
    DIAMOND_LOCATION = 'http://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz'

    input_specs: InputSpecs = InputSpecs({
        'input_path': InputSpec(File, human_name="Fasta File", short_description="Fasta Input File")
    })

    output_specs: OutputSpecs = OutputSpecs({'output_path': OutputSpec(
        Table, human_name="Blast results", short_description="This table resumes Diamond results"), })

    config_specs: ConfigSpecs = ConfigSpecs(
        {
            "input_type_value":
            StrParam(
                default_value="nuc", allowed_values=["nuc", "prot"],
                short_description="Type of alignement to perform : Prot against Prot database (i.e blastp) or Translated Nucl against prot database (i.e blastx). [Respectivly, options : prot, nuc ]. Default = nuc"),
            "evalue_value":
            FloatParam(
                default_value=0.00001, min_value=0.0, short_description="E-value : Default = 0.00001 (i.e 1e-5)"),
            "num_threads": IntParam(default_value=1, min_value=1, short_description="Number of threads"),
            "query_cover_value":
            IntParam(
                default_value=70, min_value=1, max_value=100,
                short_description="Report only alignments above the given percentage of query cover (min= 1, max= 100). [Default = 70]")})

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_blast_diamond.py"
    )

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:

        input_path: File = inputs['input_path']
        input_type_value = params["input_type_value"]
        evalue_value = params["evalue_value"]
        num_threads = params["num_threads"]
        query_cover_value = params["query_cover_value"]

        diamond_result_path = self.call_diamond(
            input_path.path, input_type_value, evalue_value, num_threads, query_cover_value)

        raw_table_columns = ['Subject ID', 'Percentage of identical matches', 'Alignment length',
                             'Number of mismatches', 'Number of gap openings', 'Start of alignment in query',
                             'End of alignment in query', 'Start of alignment in subject',
                             'End of alignment in subject', 'Expected value', 'Bit score']

        # Check if Diamond produced any hits (empty file = 0 hits)
        file_has_data = False
        with open(diamond_result_path) as f:
            for line in f:
                if not line.startswith('#') and line.strip():
                    file_has_data = True
                    break

        if not file_has_data:
            self.log_info_message("Diamond returned 0 hits. Returning empty table.")
            diamond_blast_table = Table(pd.DataFrame(columns=raw_table_columns))
        else:
            # Use TableImporter to read the raw table
            diamond_blast_table: Table = TableImporter.call(
                File(diamond_result_path),
                {'delimiter': 'tab', 'header': -1, 'file_format': 'tsv',
                    'index_column': 0, 'comment': '#'}
            )
            diamond_blast_table.set_all_column_names(raw_table_columns)

        # Return the output table
        return {
            'output_path': diamond_blast_table
        }

    def call_diamond(self, diamond_input_path: str, input_type_value: str, evalue_value: float, num_threads: int,
                     query_cover_value: int) -> str:

        file_downloader = TaskFileDownloader(brick_name=Diamond.get_brick_name(),
                                             message_dispatcher=self.message_dispatcher)

        diamond_folder_path = file_downloader.download_file_if_missing(self.DIAMOND_LOCATION, 'diamond',
                                                                       decompress_file=True)

        swissprot_db_path = file_downloader.download_file_if_missing(self.DB_LOCATION, 'swissprot',
                                                                     decompress_file=True)

        shell_proxy = ShellProxy(message_dispatcher=self.message_dispatcher)

        diamond_exe_path = os.path.join(diamond_folder_path, 'diamond')
        swissprot_exe_path = os.path.join(swissprot_db_path, 'swissprot')

        output_file = os.path.join(shell_proxy.working_dir, 'output.tsv')

        # blastp for protein input, blastx for nucleotide input
        blast_mode = 'blastp' if input_type_value == 'prot' else 'blastx'

        command = [
            diamond_exe_path,
            blast_mode,
            '--evalue', str(evalue_value),
            '--query-cover', str(query_cover_value),
            '-q', diamond_input_path,
            '-o', output_file,
            '--threads', str(num_threads),
            '-d', swissprot_exe_path,
            '--header'
        ]
        shell_proxy.run(command)

        return output_file
