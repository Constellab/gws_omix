
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (CondaShellProxy, ConfigParams, ConfigSpecs, File,
                      FileHelper, InputSpec, InputSpecs, IntParam, OutputSpec,
                      OutputSpecs, Settings, Table, TableImporter, Task,
                      TaskFileDownloader, TaskInputs, TaskOutputs,
                      task_decorator)


@task_decorator("DeepECtransformer", human_name="DeepECtransformer",
                short_description="Deep ec v2")
class DeepECtransformer(Task):

    input_specs: InputSpecs = InputSpecs({
        'input_path': InputSpec(File, human_name="Fasta File", short_description="Fasta Input File")
    })

    output_specs: OutputSpecs = OutputSpecs({'output_table': OutputSpec(
        Table, human_name="DeepEC results", short_description="This table resumes DeepEC results"), })

    config_specs: ConfigSpecs = ConfigSpecs({
        "num_threads": IntParam(default_value=4, min_value=1, short_description="Number of threads"),
    })

    SOURCE_CODE = "https://github.com/kaistsystemsbiology/DeepProZyme/archive/refs/tags/v_1_0.zip"
    DEEP_EC_FOLDER_NAME = 'DeepProZyme-v_1_0'
    OUTPUT_FILE_NAME = 'DeepECv2_result.txt'

    temp_dir: str = None

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:

        deepec_folder = self.download_deepec_repository()

        #######################

        input_file: File = inputs['input_path']
        num_threads = params["num_threads"]

        ec_result_path = self.call_deep_transformer(
            input_file.path, num_threads, deepec_folder)

        # Use TableImporter to read the raw table
        ec_table: Table = TableImporter.call(
            File(ec_result_path),
            {'delimiter': 'tab', 'header': 0,
                'file_format': 'txt', 'index_column': -1}
        )

        # rename columns
        raw_table_columns = ["sequence_ID", "prediction"]
        ec_table.set_all_column_names(raw_table_columns)

        # Return the output table
        return {
            'output_table': ec_table
        }

    def download_deepec_repository(self) -> str:
        file_downloader = TaskFileDownloader(brick_name=DeepECtransformer.get_brick_name(),
                                             message_dispatcher=self.message_dispatcher)
        return file_downloader.download_file_if_missing(self.SOURCE_CODE, 'deepECv2', decompress_file=True)

    def call_deep_transformer(self, input_path: str, num_threads: int, deepec_folder: str) -> str:

        repo_path = os.path.join(deepec_folder, self.DEEP_EC_FOLDER_NAME)
        python_file = os.path.join(
            deepec_folder, self.DEEP_EC_FOLDER_NAME, "run_deepectransformer.py")
        yml_env_file = os.path.join(
            deepec_folder, self.DEEP_EC_FOLDER_NAME, "environment.yml")

        # we set the working to the repository, because to work the command needs to be run inside the repository
        shell_proxy = CondaShellProxy(env_file_path=yml_env_file, env_name="DeepProZyme",
                                      message_dispatcher=self.message_dispatcher, working_dir=repo_path)

        # create a temp dir to store output in this
        self.temp_dir = Settings.make_temp_dir()
        result = shell_proxy.run(["python", python_file, "-i", input_path, '-o',
                                  self.temp_dir, '-g', "cpu", "-cpu", str(num_threads)])

        if result != 0:
            raise Exception('Error during deepec command')

        return os.path.join(self.temp_dir, self.OUTPUT_FILE_NAME)

    def run_after_task(self):
        if self.temp_dir:
            FileHelper.delete_dir(self.temp_dir)


# "python run_deepectransformer.py -i ./example/mdh_ecoli.fa -o ./example/results -g cpu -b 128 -cpu 2"
