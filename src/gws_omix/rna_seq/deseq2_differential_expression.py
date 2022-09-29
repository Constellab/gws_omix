# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import glob
import os

from gws_core import (BadRequestException, File, Folder, InputSpec, IntParam,
                      MetadataTable, MetadataTableImporter, OutputSpec,
                      Settings, StrParam, Table, TableImporter,
                      TableRowAnnotatorHelper, TaskInputs, TaskOutputs,
                      task_decorator)
from gws_core.config.config_types import ConfigParams, ConfigSpecs
from gws_core.io.io_spec import InputSpec, OutputSpec
from gws_core.io.io_spec_helper import InputSpecs, OutputSpecs
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
    input_specs: InputSpecs = {
        'salmon_reads_quantmerge_file':
        InputSpec(
            SalmonReadsQuantmergeOutputFile, human_name="Salmon_merged_counts",
            short_description="Salmon merged raw count files"),
        'metadata_file':
        InputSpec(
            File, human_name="Metadata_file",
            short_description="Metadata file describing samples (see https://hub.gencovery.com/bricks/gws_omix/latest/doc/use-cases/undefined )")}
    output_specs: OutputSpecs = {'DESeq2_tables': OutputSpec(
        ResourceSet, human_name="Deseq2 output files", short_description="DEseq2 output tables containing fold-change and statistics")}
    config_specs: ConfigSpecs = {
        "metadata_column": StrParam(
            human_name="Metadata column",
            optional=True,
            short_description="Column on which the differential analysis will be performed"),
        # "design_formula": StrParam(human_name="Extended experimental design formula",
        #                            short_description="For complex multifactorial analyses, add extra metadata column. ex: Size + Time + Time:Size (see: https://cran.r-project.org/doc/manuals/R-intro.html#Formulae-for-statistical-models )",
        #                            optional=True, visibility=StrParam.PROTECTED_VISIBILITY, default_value=None),
        "output_file_names": StrParam(
            default_value="DESeq2",
            human_name="Output file names",
            short_description="Choose the output file names (e.g. output_file --> output-file.XXX.txt, ...)")
    }

    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        # Get output folder
        #result_folder = Folder()
        #result_folder.path = self.working_dir
        # Create ressource set containing differnetial analysis results tables
        resource_table_set: ResourceSet = ResourceSet()
        resource_table_set.name = "Set of DESeq2 differential analysis tables"
        for path in glob.glob(os.path.join(self.working_dir, "*.txt")):
            #path = os.path.join(path, value)
            table = TableImporter.call(File(path=path), {'delimiter': 'tab', "index_column": 0})
            basename = os.path.basename(path)
            table.name = basename
            # if params["design_formula"]:
            #     table.tags = {"Design": params["design_formula"]}
            # else:
            #     table.tags = {"Design": params["metadata_column"]}
            resource_table_set.add_resource(table)
        return {
            'DESeq2_tables': resource_table_set
        }

    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        salmon_merged_matrix = inputs["salmon_reads_quantmerge_file"]
        metadata = inputs["metadata_file"]
        metadata_col = params["metadata_column"]
        output_file_id = params["output_file_names"]
        #design_formula = params["design_formula"]
        script_file_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [
            " Rscript --vanilla ",
            #os.path.join(script_file_dir, "./R/Deseq2_script.one_parameter.R"),
            os.path.join(script_file_dir, "./R/deseq2.venn_file.R"),
            salmon_merged_matrix.path,
            metadata.path,
            metadata_col,
            output_file_id
        ]

        return cmd

        # if params["design_formula"]:
        #     if params["metadata_column"]:
        #         raise BadRequestException("Only metadata column or design formula is required")
        #     design_formula = params["design_formula"]
        #     cmd = [
        #         " Rscript --vanilla ",
        #         os.path.join(script_file_dir, "./R/Deseq2_script.multi_parameter_interaction.R"),
        #         salmon_merged_matrix.path,
        #         metadata.path,
        #         "\' " + design_formula + " \' ",
        #         output_file_id
        #     ]
        # elif params["metadata_column"]:
        #     metadata_col = params["metadata_column"]
        #     cmd = [
        #         " Rscript --vanilla ",
        #         os.path.join(script_file_dir, "./R/Deseq2_script.one_parameter.R"),
        #         salmon_merged_matrix.path,
        #         metadata.path,
        #         metadata_col,
        #         output_file_id
        #     ]
        # else:
        #     raise BadRequestException("Metadata column or design formula is required")

        # return cmd
