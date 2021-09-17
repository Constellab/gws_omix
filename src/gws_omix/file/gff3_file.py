# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from gws_core import File, resource_decorator

@resource_decorator("GFF3File",
                    human_name="GFF3File",
                    short_description="GFF3 File")
class GFF3File(File):
    ''' GFF3 file class'''
