# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from gws_core import File, resource_decorator

@resource_decorator("salmon_reads_quantmerge_output_file",
                    human_name="salmon_reads_quantmerge_output_file",
                    short_description="salmon_reads_quantmerge_output_file")
class SalmonReadsQuantmergeOutputFile(File):
    """salmon_Reads_quantmerge_output_file class"""
