# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from gws_core import File, resource_decorator

@resource_decorator("FastqFile",
                    human_name="FastqFile",
                    short_description="Fastq File")


class FastqFile(File):
    ''' Fastq file class'''
