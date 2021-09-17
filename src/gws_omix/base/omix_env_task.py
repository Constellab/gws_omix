# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
from gws_core import task_decorator, CondaEnvShell

@task_decorator("BaseOmixEnvTask")
class BaseOmixEnvTask(CondaEnvShell):
    __cdir = os.path.abspath(os.path.dirname(__file__))
    env_file_path = os.path.join(__cdir, "omix_env.yml")
