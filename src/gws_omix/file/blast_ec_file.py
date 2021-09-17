# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from gws_core import File, resource_decorator

@resource_decorator("BlastECFile",
                    human_name="BlastECFile",
                    short_description="BlastEC File")
class BlastECFile(File):
    ''' BlastEC file class'''