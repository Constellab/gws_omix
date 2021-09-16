# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from gws_core import File, resource_decorator

@resource_decorator("BlastEcFile",
                    human_name="BlastEcFile",
                    short_description="BlastEc File")

class BlastEcFile(File):
    ''' BlastEc file class'''