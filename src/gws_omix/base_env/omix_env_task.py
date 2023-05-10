# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import CondaEnvShell, CondaShellProxy, task_decorator
from gws_core.core.classes.observer.message_dispatcher import MessageDispatcher


def create_omix_conda_env(working_dir: str = None) -> CondaEnvShell:
    return CondaShellProxy(BaseOmixEnvTask.unique_env_name,
                           BaseOmixEnvTask.env_file_path, working_dir)


@task_decorator("BaseOmixEnvTask", hide=True)
class BaseOmixEnvTask(CondaEnvShell):
    unique_env_name = "BaseOmixEnvTask"
    env_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "omix_env.yml"
    )


class BaseOmixEnvHelper():
    unique_env_name = "BaseOmixEnvTask"
    env_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "omix_env.yml"
    )

    @classmethod
    def create_proxy(cls, message_dispatcher: MessageDispatcher = None) -> CondaShellProxy:
        return CondaShellProxy(cls.unique_env_name, cls.env_file_path, message_dispatcher=message_dispatcher)
