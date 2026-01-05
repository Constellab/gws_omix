# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (
    ConfigParams,
    Folder,
    InputSpec,
    InputSpecs,
    OutputSpec,
    OutputSpecs,
    Task,
    TaskInputs,
    TaskOutputs,
    task_decorator,
)

from .multiqc_env import MultiQcShellProxyHelper


@task_decorator(
    "MultiQC",
    human_name="MultiQC",
    short_description="Aggregate FastQC reports into a single MultiQC report.",
)
class MultiQc(Task):
    """Run MultiQC on a folder containing FastQC outputs and return a combined report."""

    input_specs: InputSpecs = InputSpecs({
        "fastqc_reports_folder": InputSpec(
            Folder,
            human_name="FastQC reports folder",
            short_description="Folder containing FastQC .html / .zip outputs",
        )
    })

    output_specs: OutputSpecs = OutputSpecs({
        "output": OutputSpec(
            Folder,
            human_name="Combined quality report",
            short_description="MultiQC combined report (HTML + data)",
        )
    })


    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """Execute MultiQC and return the output folder."""

        # Retrieve input folder path
        src = inputs["fastqc_reports_folder"].path  # type: str

        if not os.path.exists(src):
            raise FileNotFoundError(f"Input folder not found: {src}")

        # Create working directory
        shell_proxy = MultiQcShellProxyHelper.create_proxy(self.message_dispatcher)
        result_dir = os.path.join(shell_proxy.working_dir, "multiqc_result")
        os.makedirs(result_dir, exist_ok=True)

        # Build command: pass the folder to MultiQC so it autodetects FastQC files
        cmd = f"multiqc {src} -o {result_dir} -n multiqc_combined.html"

        retcode = shell_proxy.run(cmd, shell_mode=True)
        if retcode != 0:
            raise RuntimeError("MultiQC finished with a non-zero exit status (see logs).")

        return {"output": Folder(result_dir)}
