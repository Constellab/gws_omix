# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
from gws_core import task_decorator, CondaEnvShell

@task_decorator("TrimFqEnvTask", hide=True)
class TrimFqEnvTask(CondaEnvShell):
    unique_env_name = "TrimFqEnvTask"
    env_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)), 
        "trimfq_env.yml"
    )
