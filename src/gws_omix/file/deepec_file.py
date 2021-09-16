# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from gws_core import File, resource_decorator

@resource_decorator("DeepEcFile",
                    human_name="DeepEcFile",
                    short_description="DeepEc File")



class DeepEcFile(File):
    ''' DeepEc file class'''