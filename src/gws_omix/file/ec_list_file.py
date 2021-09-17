# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from gws_core import File, resource_decorator

@resource_decorator("ECListFile",
                    human_name="ECListFile",
                    short_description="File of EC list")
class ECListFile(File):
    ''' EC list file class'''