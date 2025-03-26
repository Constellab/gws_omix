# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, Folder, IntParam, Task, InputSpecs, OutputSpecs,
                      TaskInputs, TaskOutputs, task_decorator, InputSpec, OutputSpec, ConfigSpecs, File)


from .star_index_env import StarIndexShellProxyHelper


@task_decorator("StarIndex", human_name="Building a genome index",
                short_description="Indexation of genome using STAR")
class StarIndex(Task):
    """
    Generating genome indexes files.In this step user supplied the reference genome sequences (FASTA files) and annotations (GTF file), from which STAR generate genome indexes that are utilized in mapping step such STAR solo to generate count matrix for scRNA seq data.

    """

    input_specs: InputSpecs = InputSpecs({
        'genome_file': InputSpec(File, short_description="Fasta file", human_name="Fasta file"),
        'annotation_file': InputSpec(File, short_description="Annotation file", human_name="GTF annotation file")
    })

    output_specs: OutputSpecs = OutputSpecs({
        'output': OutputSpec(Folder, human_name="indexed genome", short_description="indexed genome")
    })

    config_specs: ConfigSpecs = ConfigSpecs({
        "threads": IntParam(default_value=2, min_value=2, short_description="Number of threads ")
    })

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # retrive the input table
        fasta_input_file: File = inputs['genome_file']
        GTF_input_file: File = inputs['annotation_file']
        thrd = params["threads"]

        # retrieve the factor param value
        shell_proxy = StarIndexShellProxyHelper.create_proxy(
            self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'result')
        os.makedirs(result_path)

        star_cmd = f'STAR --runMode genomeGenerate --runThreadN {thrd} --genomeChrBinNbits 12 --limitGenomeGenerateRAM 60000000000 --genomeDir {result_path} --genomeFastaFiles {fasta_input_file.path} --sjdbGTFfile {GTF_input_file.path} --genomeSAsparseD 3'
        res = shell_proxy.run(star_cmd, shell_mode=True)

        if res != 0:
            raise Exception("One error occured when formating output files")

        folder = Folder(result_path)

        # return the output table
        return {'output': folder}
