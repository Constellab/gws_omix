# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (CondaEnvShell, CondaShellProxy, MessageDispatcher,
                      task_decorator)


@task_decorator("PyGSEAEnvTask", hide=True)
class PyGSEAEnvTask(CondaEnvShell):
    unique_env_name = "PyGSEAEnvTask"
    env_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "pygsea_env.yml"
    )


class PyGSEAShellProxyHelper():
    ENV_DIR_NAME = "PyGSEAShellProxy"
    ENV_FILE_PATH = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "pygsea_env.yml"
    )

    # @classmethod
    # def create_proxy(cls):
    #     return CondaShellProxy(cls.ENV_DIR_NAME, cls.ENV_FILE_PATH)
    @classmethod
    def create_proxy(cls, message_dispatcher: MessageDispatcher = None):
        return CondaShellProxy(cls.ENV_DIR_NAME, cls.ENV_FILE_PATH, message_dispatcher=message_dispatcher)
