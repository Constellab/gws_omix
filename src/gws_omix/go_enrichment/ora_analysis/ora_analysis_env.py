# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import MessageDispatcher, PipShellProxy


class OraAnalysisShellProxyHelper:
    ENV_DIR_NAME = "OraAnalysisShellProxy2"
    ENV_FILE_PATH = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "ora_analysis_env.txt"
    )

    @classmethod
    def create_proxy(cls, message_dispatcher: MessageDispatcher = None):
        return PipShellProxy(cls.ENV_FILE_PATH, cls.ENV_DIR_NAME, message_dispatcher=message_dispatcher)
