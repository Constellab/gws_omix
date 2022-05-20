# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from gws_core import Folder, resource_decorator


@resource_decorator("StarIndexFolder", human_name="STAR Index Folder",
                    short_description="STAR Index Result Folder")
class StarIndexFolder(Folder):
    ''' StarIndexFolder folder class'''
