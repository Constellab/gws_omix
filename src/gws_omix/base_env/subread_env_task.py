# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import CondaShellProxy, MessageDispatcher


class SubreadShellProxyHelper():
    ENV_DIR_NAME = "SubreadShellProxy"
    ENV_FILE_PATH = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "./subread_env.yml"
    )

    @classmethod
    def create_proxy(cls, message_dispatcher: MessageDispatcher = None):
        return CondaShellProxy(cls.ENV_DIR_NAME, cls.ENV_FILE_PATH, message_dispatcher=message_dispatcher)
