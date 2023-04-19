# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import MessageDispatcher, PipShellProxy

__cdir__ = os.path.dirname(os.path.realpath(__file__))


class PygseaPipShellProxyHelper():
    ENV_DIR_NAME = "PygseaPipShellProxyHelper"
    ENV_FILE_PATH = os.path.join(
        __cdir__,
        "pygsea_pip_env.txt"
    )

    @classmethod
    def create_proxy(cls, message_dispatcher: MessageDispatcher = None):
        return PipShellProxy(cls.ENV_DIR_NAME, cls.ENV_FILE_PATH,
                             message_dispatcher=message_dispatcher)

# @task_decorator("PygseaPipEnvTask")
# class PygseaPipEnvTask(PipEnvTask):
#     input_specs = {}
#     output_specs = {'file': OutputSpec(File)}

#     unique_env_name = "PygseaPipEnvTask"
#     env_file_path = os.path.join(
#         __cdir__,
#         "pygsea_pip_env.yml"
#     )

#     async def run_with_proxy(self, params: ConfigParams, inputs: TaskInputs, shell_proxy: PipShellProxy) -> TaskOutputs:
#         command = ["python", os.path.join(__cdir__, "penv", "jwt_encode.py"), ">", "out.txt"]
#         shell_proxy.run(command, shell_mode=True)

#         # retrieve the result
#         file = File(path=os.path.join(self.working_dir, "out.txt"))
#         return {"file": file}

# class PygseaPipEnvProxyHelper():
#     ENV_DIR_NAME = "PygseaPipEnvlProxy"
#     ENV_FILE_PATH = os.path.join(
#         os.path.abspath(os.path.dirname(__file__)),
#         "pygsea_pip_env.yml"
#     )

#     @classmethod
#     def create_proxy(cls, message_dispatcher: MessageDispatcher = None):
#         return PipShellProxy(cls.ENV_DIR_NAME, cls.ENV_FILE_PATH, message_dispatcher=message_dispatcher)
