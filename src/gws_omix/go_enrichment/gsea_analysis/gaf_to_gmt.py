# gaf_to_gmt.py
# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import shlex
from pathlib import Path
from typing import Final

from gws_core import (
    ConfigParams,
    ConfigSpecs,
    File,
    InputSpec,
    InputSpecs,
    OutputSpec,
    OutputSpecs,
    ResourceSet,
    ShellProxy,
    StrParam,
    Task,
    TaskInputs,
    TaskOutputs,
    task_decorator,
)

from .gsea_env import GseaShellProxyHelper


@task_decorator(
    "GafToGmt",
    human_name="GAF → GMT (GO_ALL + BP/MF/CC; withIEA/noIEA × id/symbol)",
    short_description="Builds GMT files from a GO GAF: GO_ALL (4 files) + split BP/MF/CC.",hide=True
)
class GafToGmt(Task):
    """
    Convert GO annotations in GAF format into GMT gene set files.

    Outputs:
      1) gmt_files (existing): GO_ALL × withIEA/noIEA × id/symbol (4 GMT)
      2) gmt_withIEA_by_go: BP/MF/CC × id/symbol (6 GMT, with IEA)
      3) gmt_noIEA_by_go:   BP/MF/CC × id/symbol (6 GMT, excluding IEA)

    Notes:
      - 'NOT' annotations are always excluded.
      - No term-size filtering is applied here (keep GMT complete).
    """

    input_specs: Final[InputSpecs] = InputSpecs({
        "gaf_file": InputSpec(
            File,
            human_name="GO Annotations (GAF/GAF.GZ)",
            short_description="Example: sgd.gaf.gz (yeast), fb.gaf.gz (drosophila).",
        ),
    })

    output_specs: Final[OutputSpecs] = OutputSpecs({
        "gmt_files": OutputSpec(ResourceSet, human_name="GMT (GO_ALL)"),
        "gmt_withIEA_by_go": OutputSpec(ResourceSet, human_name="GMT (withIEA) split by GO aspect"),
        "gmt_noIEA_by_go": OutputSpec(ResourceSet, human_name="GMT (noIEA) split by GO aspect"),
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "prefix": StrParam(default_value="organism", short_description="Prefix for output GMT filenames."),
    })

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        gaf: File = inputs["gaf_file"]
        prefix = str(params.get("prefix", "organism")).strip() or "organism"

        shell_proxy: ShellProxy = GseaShellProxyHelper.create_proxy(self.message_dispatcher)
        work = Path(shell_proxy.working_dir)
        outdir = work / "gmt_out"
        outdir.mkdir(parents=True, exist_ok=True)

        script_dir = Path(os.path.dirname(os.path.realpath(__file__)))
        script_path = script_dir / "_gaf_to_gmt.py"

        cmd = [
            "python3",
            str(script_path),
            "--gaf", gaf.path,
            "--outdir", str(outdir),
            "--prefix", prefix,
        ]
        cmd_str = " ".join(shlex.quote(str(x)) for x in cmd)

        ret = shell_proxy.run(cmd_str, shell_mode=True)
        if ret != 0:
            return {
                "gmt_files": ResourceSet(),
                "gmt_withIEA_by_go": ResourceSet(),
                "gmt_noIEA_by_go": ResourceSet(),
            }

        # --- RS1: existing GO_ALL (4) ---
        rs_all = ResourceSet()
        rs_all.name = "GMT (GO_ALL) — withIEA/noIEA × id/symbol"

        expected_all = {
            "withIEA_id": outdir / f"{prefix}_GO_ALL_withIEA_id.gmt",
            "withIEA_symbol": outdir / f"{prefix}_GO_ALL_withIEA_symbol.gmt",
            "noIEA_id": outdir / f"{prefix}_GO_ALL_noIEA_id.gmt",
            "noIEA_symbol": outdir / f"{prefix}_GO_ALL_noIEA_symbol.gmt",
        }
        for key, fp in expected_all.items():
            if fp.exists() and fp.stat().st_size > 0:
                rs_all.add_resource(File(path=str(fp)), key)

        # --- RS2: withIEA split by aspect (BP/MF/CC × id/symbol) ---
        rs_with = ResourceSet()
        rs_with.name = "GMT (withIEA) — split by GO aspect (BP/MF/CC) × id/symbol"

        expected_with = {}
        for aspect in ("BP", "MF", "CC"):
            expected_with[f"{aspect}_id"] = outdir / f"{prefix}_GO_{aspect}_withIEA_id.gmt"
            expected_with[f"{aspect}_symbol"] = outdir / f"{prefix}_GO_{aspect}_withIEA_symbol.gmt"

        for key, fp in expected_with.items():
            if fp.exists() and fp.stat().st_size > 0:
                rs_with.add_resource(File(path=str(fp)), key)

        # --- RS3: noIEA split by aspect (BP/MF/CC × id/symbol) ---
        rs_no = ResourceSet()
        rs_no.name = "GMT (noIEA) — split by GO aspect (BP/MF/CC) × id/symbol"

        expected_no = {}
        for aspect in ("BP", "MF", "CC"):
            expected_no[f"{aspect}_id"] = outdir / f"{prefix}_GO_{aspect}_noIEA_id.gmt"
            expected_no[f"{aspect}_symbol"] = outdir / f"{prefix}_GO_{aspect}_noIEA_symbol.gmt"

        for key, fp in expected_no.items():
            if fp.exists() and fp.stat().st_size > 0:
                rs_no.add_resource(File(path=str(fp)), key)

        return {
            "gmt_files": rs_all,
            "gmt_withIEA_by_go": rs_with,
            "gmt_noIEA_by_go": rs_no,
        }