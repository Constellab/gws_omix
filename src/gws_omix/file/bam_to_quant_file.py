# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from gws_core import File, resource_decorator

@resource_decorator("BAMToQuantFile",
                    human_name="BAMToQuantFile",
                    short_description="BAMToQuantFile")
class BAMToQuantFile(File):
    """BAMToQuantFile file class"""
