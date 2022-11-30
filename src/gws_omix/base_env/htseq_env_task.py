# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
from gws_core import task_decorator, CondaEnvShell, CondaShellProxy


@task_decorator("HTSeqEnvTask", hide=True)
class HTSeqEnvTask(CondaEnvShell):
    unique_env_name = "HTSeqEnvTask"
    env_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "htseq_env.yml"
    )


class HTSeqShellProxyHelper():
    ENV_DIR_NAME = "HTSeqShellProxy"
    ENV_FILE_PATH = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "htseq_env.yml"
    )

    @classmethod
    def create_proxy(cls):
        return CondaShellProxy(cls.ENV_DIR_NAME, cls.ENV_FILE_PATH)
