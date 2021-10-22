# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
from gws_core import task_decorator, CondaEnvShell

@task_decorator("TrimFqEnvTask")
class TrimFqEnvTask(CondaEnvShell):
    __cdir = os.path.abspath(os.path.dirname(__file__))
    env_file_path = os.path.join(__cdir, "trimfq_env.yml")
