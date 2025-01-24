# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import MessageDispatcher, PipShellProxy


class HTSeqShellProxyHelper():

    ENV_DIR_NAME = "HTSeqShellProxyPip"
    ENV_FILE_PATH = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "htseq_env_2.txt"
    )

    @classmethod
    def create_proxy(cls, message_dispatcher: MessageDispatcher = None):
        return PipShellProxy(env_file_path=cls.ENV_FILE_PATH, env_name=cls.ENV_DIR_NAME,
                             message_dispatcher=message_dispatcher)
