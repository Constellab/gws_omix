# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import MessageDispatcher, PipShellProxy


class PygseaPipShellProxyHelper():
    ENV_DIR_NAME = "PygseaPipShellProxyHelper"
    ENV_FILE_PATH = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "pygsea_pip_env.txt"
    )

    @classmethod
    def create_proxy(cls, message_dispatcher: MessageDispatcher = None):
        return PipShellProxy(cls.ENV_DIR_NAME, cls.ENV_FILE_PATH,
                             message_dispatcher=message_dispatcher)
