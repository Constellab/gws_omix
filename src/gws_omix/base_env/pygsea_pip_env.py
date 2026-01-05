

import os

from gws_core import MessageDispatcher, PipShellProxy


class PygseaPipShellProxyHelper:
    ENV_DIR_NAME = "PygseaPipShellProxyHelper"
    ENV_FILE_PATH = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "pygsea_pip_env.txt"
    )

    @classmethod
    def create_proxy(cls, message_dispatcher: MessageDispatcher = None):
        return PipShellProxy(env_file_path=cls.ENV_FILE_PATH, env_name=cls.ENV_DIR_NAME,
                             message_dispatcher=message_dispatcher)
