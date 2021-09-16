# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from gws_core import File, resource_decorator

@resource_decorator("EcListFile",
                    human_name="EcListFile",
                    short_description="EcList File")



class EcListFile(File):
    ''' EcList file class'''