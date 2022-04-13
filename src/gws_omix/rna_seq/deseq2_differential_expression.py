# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import glob
import os

from gws_core import (ConfigParams, File, IntParam, MetadataTable,
                      MetadataTableImporter, Settings, StrParam, Table,
                      TableImporter, TableRowAnnotatorHelper, TaskInputs,
                      TaskOutputs, task_decorator)
from gws_core.resource.resource_set import ResourceSet

from ..base_env.r_env_task import BaseREnvTask
from ..file.salmon_reads_quantmerge_output_file import \
    SalmonReadsQuantmergeOutputFile


@task_decorator("DESeq2DifferentialAnalysis", human_name="DESeq2 pairwise differential analysis",
                short_description="Compute differential analysis using DESeq2 R package (pairwise comparison)")
class DESeq2DifferentialAnalysis(BaseREnvTask):
    """
    DESeq2DifferentialAnalysis class.
    """

    input_specs = {
        'salmon_reads_quantmerge_file': SalmonReadsQuantmergeOutputFile,
        'metadata_file': File
    }
    output_specs = {
        'DESeq2_tables': ResourceSet
    }
    config_specs = {
        # "taxonomic_level":
        # IntParam(
        #     min_value=1, human_name="Taxonomic level",
        #     short_description="Taxonomic level id: 1_Kingdom, 2_Phylum, 3_Class, 4_Order,5_Family, 6_Genus, 7_Species"),
        "metadata_column": StrParam(
            human_name="Metadata column",
            short_description="Column on which the differential analysis will be performed"),
        "output_file_names": StrParam(
            default_value="DESeq2_output_file",
            human_name="Output file names",
            short_description="Choose the output file names (e.g. output_file --> output-file.XXX.txt, ...)")
        # "threads": IntParam(default_value=2, min_value=2, short_description="Number of threads")
    }

    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:

        #  Importing Metadata table
        path = inputs["metadata_file"]
        metadata_table = MetadataTableImporter.call(File(path=path), {'delimiter': 'tab'})

        # Create ressource set containing differnetial analysis results tables
        resource_table_set: ResourceSet = ResourceSet()
        resource_table_set.name = "Set of DESeq2 differential analysis tables"
        # deseq2_output_files = os.path.join(self.working_dir, "*.txt")

        # deseq2_output_file_path = []
        for output_file_path in glob.glob(os.path.join(self.working_dir, "*.txt")):
            # deseq2_output_file_path.append(output_file_path)

            # for output_file_path in self.deseq2_output_file_paths.items():
            table = TableImporter.call(File(path=output_file_path), {'delimiter': 'tab', "index_column": 0})
            table_annotated = TableRowAnnotatorHelper.annotate(table, metadata_table)
            basename = os.path.basename(output_file_path)
            table_annotated.name = basename
            resource_table_set.add_resource(table_annotated)

        return {
            'DESeq2_tables': resource_table_set
        }

    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        salmon_merged_matrix = inputs["salmon_reads_quantmerge_file"]
        metadata = inputs["metadata_file"]
        metadata_col = params["metadata_column"]
        output_file_id = params["output_file_names"]
        # thrds = params["threads"]
        script_file_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [
            " Rscript --vanilla ",
            os.path.join(script_file_dir, "./R/DESeq2_salmon_differential_expression.R"),
            salmon_merged_matrix.path,
            metadata.path,
            metadata_col,
            output_file_id
        ]
        return cmd

    # def _get_output_file_path(self):
    #     return os.path.join(
    #         self.working_dir,
    #         "diversity"
    #     )
