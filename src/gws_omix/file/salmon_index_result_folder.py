# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from gws_core import Folder, resource_decorator


@resource_decorator("SalmonIndexResultFolder", human_name="SalmonIndexResultFolder",
                    short_description="Salmon Index Result Folder")
class SalmonIndexResultFolder(Folder):
    ''' SalmonIndexResultFolder file class'''
