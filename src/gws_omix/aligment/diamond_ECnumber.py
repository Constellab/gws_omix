import os
from gws_core import (ConfigParams, ConfigSpecs, File, ShellProxy,
                      FloatParam, InputSpec, InputSpecs, IntParam, OutputSpec, TableImporter,
                      OutputSpecs, StrParam, Task, TaskFileDownloader, Table,
                      TaskInputs, TaskOutputs, task_decorator, StrParam)

from typing import Dict, List

@task_decorator("Ec_number", human_name="Ec_number",
                short_description="xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx.")

class ECnumber(Task):
    """
    BlastEC class.

    Represents a process that wraps NCBI blast program. This version
    allows getting EC numbers for digital twins reconstruction.
    """

    EC_LOCATION = 'https://storage.gra.cloud.ovh.net/v1/AUTH_a0286631d7b24afba3f3cdebed2992aa/opendata/omix/db/fungi_EC_number_database_uniprot.tsv'
    DB_FILENAME = 'fungi_EC_number_database_uniprot.tsv'

    input_specs: InputSpecs = InputSpecs({
        'input_table': InputSpec(Table, human_name="diamond result", short_description="Diamond result")
    })

    output_specs: OutputSpecs = OutputSpecs({
        'output_table': OutputSpec(Table, human_name="xxxxxxxxxxxxx", short_description="xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
    })

    config_specs: ConfigSpecs = {}

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:

        input_table: Table = inputs['input_table']
        ec_blast_table = self.call_ecnumber(input_table)

        return {
            'output_table': ec_blast_table
        }

    def call_ecnumber(self, input_table: Table) -> Table:
        file_downloader = TaskFileDownloader(
            brick_name=ECnumber.get_brick_name(),
            message_dispatcher=self.message_dispatcher
        )

        ec_database_file = file_downloader.download_file_if_missing(
            self.EC_LOCATION, filename=self.DB_FILENAME, decompress_file=False
        )

        intable = input_table.get_column_data("Subject ID")
        processed_data = []
        for value in intable:
            # Extract information from each value in the "Subject ID" column and keep only the integer section
            data = value.split('.')[0]
            processed_data.append(data)

        # read the databse file
        mapping_dict : Dict[str, List[str]] = {}
        with open(ec_database_file, 'r') as mapping_file:
            header = mapping_file.readline().strip().split('\t')
            entry_index = header.index('Entry')
            ec_number_index = header.index('EC number')

            for line in mapping_file:
                columns = line.strip().split('\t')
                entry = columns[entry_index]
                ec_numbers = columns[ec_number_index].split(';')
                mapping_dict[entry] = ec_numbers


        # compare database file to diamond table
        processed_ec_numbers = set()
        # Create a list to store the processed EC numbers

        for entry in processed_data:
            if entry in mapping_dict and mapping_dict[entry]:
                for ec_number in mapping_dict[entry]:
                    ec_number_stripped = ec_number.strip()
                    if ec_number_stripped and ec_number_stripped not in processed_ec_numbers:
                        processed_ec_numbers.add(ec_number_stripped)

        # Create a Table with the processed EC numbers
        result = Table(list(processed_ec_numbers),column_names= ["EC number"])
        return result