# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import re
import pandas as pd
from gws_core import (CondaShellProxy, ConfigParams, ConfigSpecs, File, ShellProxy,
                      FloatParam, InputSpec, InputSpecs, IntParam, OutputSpec, TableImporter,
                      OutputSpecs, StrParam, Task, TaskFileDownloader, Table,
                      TaskInputs, TaskOutputs, task_decorator, StrParam, TaskFileDownloader)

from gws_omix.base_env.diamond_env_task import DiamondShellProxyHelper


@task_decorator("Diamond", human_name="diamond",
                short_description="xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx.")
class Diamond(Task):
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

    DB_LOCATION = 'https://storage.gra.cloud.ovh.net/v1/AUTH_a0286631d7b24afba3f3cdebed2992aa/opendata/omix/db/swissprot_database.zip'
    DIAMOND_LOCATION = 'http://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz'

    input_specs: InputSpecs = InputSpecs({
        'input_path': InputSpec(File, human_name="Fasta File", short_description="Fasta Input File")
    })

    output_specs: OutputSpecs = OutputSpecs({'output_path': OutputSpec(
        Table, human_name="xxxxxxxxxxxxx", short_description="xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"), })

    config_specs: ConfigSpecs = {
        "input_type_value": StrParam(default_value="nuc", allowed_values=["nuc", "prot"], short_description="Type of alignement to perform : Prot against Prot database (i.e blastp) or Translated Nucl against prot database (i.e blastx). [Respectivly, options : PP, TNP ]. Default = PP"),
        "evalue_value": FloatParam(default_value=0.00001, min_value=0.0, short_description="E-value : Default = 0.00001 (i.e 1e-5)"),
        "num_threads": IntParam(default_value=1, min_value=1, short_description="Number of threads"),
        "query_cover_value": IntParam(default_value=90, min_value=1, max_value=100, short_description="query_cover_value"),
        "approx_id_value": IntParam(default_value=70, min_value=1, max_value=100, short_description="Coverage (see blast option -qcov_hsp_perc) minimum percentage threshold to exclude results [Default = 70]"),
        "min_score_value": IntParam(default_value=500, min_value=1, max_value=1000, short_description=""),


    }
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
        approx_id_value = params["approx_id_value"]
        min_score_value = params["min_score_value"]

        diamond_input_path: str = None
        if input_type_value == "nuc":
            diamond_input_path = self.translate(input_path.path)
        else:
            diamond_input_path = input_path.path

        diamond_result_path = self.call_diamond(diamond_input_path, evalue_value, num_threads, query_cover_value,
                                                approx_id_value, min_score_value)


        # Use TableImporter to read the raw table
        diamond_blast_table : Table = TableImporter.call(
            File(diamond_result_path),
            {'delimiter': 'tab', 'header': -1, 'file_format': 'tsv', 'index_column': 0, 'comment': '#'}
        )
        raw_table_columns = ['Subject ID', 'Percentage of identical matches', 'Alignment length',
                            'Number of mismatches', 'Number of gap openings', 'Start of alignment in query',
                            'End of alignment in query', 'Start of alignment in subject', 'End of alignment in subject',
                            'Expected value', 'Bit score']

        diamond_blast_table.set_all_column_names(raw_table_columns)

        # Return the output table
        return {
            'output_path': diamond_blast_table
        }

    def translate(self, input_file_path: str) -> str:
        # Execute the command
        shell_proxy: ShellProxy = DiamondShellProxyHelper.create_proxy(self.message_dispatcher)

        fasta_output_path = os.path.join(shell_proxy.working_dir, 'output.fasta')

        # call python file
        cmd = f"python3 {self.python_file_path} {input_file_path} {fasta_output_path}"
        shell_proxy.run(cmd, shell_mode=True)

        return fasta_output_path


    def call_diamond(self, diamond_input_path: str, evalue_value: float, num_threads: int,
                     query_cover_value: int, approx_id_value: int, min_score_value: int) -> str:

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

        command = [
            diamond_exe_path,
            'blastp',
            '--min-score', str(min_score_value),
            '--evalue', str(evalue_value),
            '--query-cover', str(query_cover_value),
            '--approx-id', str(approx_id_value),
            '-q', diamond_input_path,
            '-o', output_file,
            '--threads', str(num_threads),
            '-d', swissprot_exe_path,  # Add 'swissprot' to the prefix
            '--header'
        ]
        shell_proxy.run(command)

        return output_file