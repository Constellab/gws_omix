# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import CondaShellProxy, MessageDispatcher


class DiamondShellProxyHelper():
    unique_env_name = "DiamondShellProxy"
    env_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "Diamond_env.yml"
    )

    @classmethod
    def create_proxy(cls, message_dispatcher: MessageDispatcher = None) -> CondaShellProxy:
        return CondaShellProxy(cls.unique_env_name, cls.env_file_path, message_dispatcher=message_dispatcher)
