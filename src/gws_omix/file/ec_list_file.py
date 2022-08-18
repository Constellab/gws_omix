# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

from gws_core import (BadRequestException, File, IntParam,
                      TextView, resource_decorator, view)



@resource_decorator("ECListFile", human_name="ECListFile", short_description="File of EC list")
class ECListFile(File):
    ''' EC list file class'''

    @view(view_type=TextView, human_name="TextView", short_description="View of the EC list file ")
    def view_as_raw_text(self, **kwargs) -> dict:
        text = self.read()
        return TextView(data=text, **kwargs)
