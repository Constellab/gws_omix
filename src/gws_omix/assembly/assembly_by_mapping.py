# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (InputSpec, IntParam, OutputSpec,
                      StrParam, TaskInputs, TaskOutputs, task_decorator, ConfigParams, ConfigSpecs, 
                      InputSpec, OutputSpec, InputSpecs, OutputSpecs, ResourceSet)

from ..base_env.omix_env_task import BaseOmixEnvTask
from ..file.bam_file import BAMFile
from ..file.fasta_file import FastaFile


@task_decorator("GenomeByMappingAssembly", human_name="Genome by mapping assembly pipeline",
                short_description="Use sequencing data mapped on reference genome to assembled this unassembled genome")
class GenomeByMappingAssembly(BaseOmixEnvTask):
    """
    GenomeByMappingAssembly class.
    """

    ASSEMBLY_FILES_PATH = {
        "chromosome fasta file": "mapping_assembly.fna",
        "vcf variant file": "variants.vcf",
    }

    # Taxo stacked CSS barplot

    input_specs: InputSpecs = {
        'genome_fasta': InputSpec(FastaFile, human_name="", short_description=""),
        'mapping_file': InputSpec(BAMFile, human_name="", short_description="")
    }
    output_specs: OutputSpecs = {
        'assembly_folder': OutputSpec(ResourceSet, human_name="", short_description="")
    }
    config_specs:  ConfigSpecs = {
        "threads": IntParam(default_value=2, min_value=2, short_description="Number of threads"),
        "mapping_quality_filter":
        IntParam(
            default_value=11, min_value=0, max_value=60,
            short_description="Threshold to filter low quality mapped reads (0: No filter, 60: perfectly mapped)"),
        "chromosome_name": StrParam(short_description="Chromosome name to be assembled")}

    def gather_outputs(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:

       # Create ressource set containing assembly output files (fasta and variant calling file)

        assembly_resource_table_set: ResourceSet = ResourceSet()
        assembly_resource_table_set.name = "Containing chromosome assembly output files (fasta and variant calling files)"
        for key, value in self.ASSEMBLY_FILES_PATH.items():
            file = os.path.join(self.working_dir, value)
            file.name = key
            assembly_resource_table_set.add_resource(file)

        return {
            'assembly_folder': assembly_resource_table_set
        }

    def build_command(self, params: ConfigParams, inputs: TaskInputs) -> list:
        fasta_file = inputs["genome_fasta"]
        bam_file = inputs["mapping_file"]
        chr_name = params["chromosome_name"]
        threads = params["threads"]
        mapq_threshold = params["mapping_quality_filter"]

        script_file_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [
            " bash ",
            os.path.join(script_file_dir, "./sh/assembly_mapping.sh"),
            fasta_file.path,
            chr_name,
            bam_file.path,
            mapq_threshold,
            threads
        ]
        return cmd
