# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import CondaEnvShell, task_decorator


@task_decorator("InterProScanEnvTask", hide=True)
class InterProScanEnvTask(CondaEnvShell):
    unique_env_name = "InterProScanEnvTask"
    env_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "interproscan_env.yml"
    )
