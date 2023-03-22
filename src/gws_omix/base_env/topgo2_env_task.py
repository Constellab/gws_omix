# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (CondaEnvShell, CondaShellProxy, MessageDispatcher,
                      task_decorator)


@task_decorator("TopGO2EnvTask", hide=True)
class TopGO2EnvTask(CondaEnvShell):
    unique_env_name = "TopGO2EnvTask"
    env_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "topgo2_env.yml"
    )


class TopGO2ShellProxyHelper():
    ENV_DIR_NAME = "TopGO2ShellProxy"
    ENV_FILE_PATH = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "topgo2_env.yml"
    )

    @classmethod
    def create_proxy(cls, message_dispatcher: MessageDispatcher = None):
        return CondaShellProxy(cls.ENV_DIR_NAME, cls.ENV_FILE_PATH, message_dispatcher=message_dispatcher)
