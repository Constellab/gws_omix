# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com
'''
import os
import re

from gws_core import (
    ConfigParams, Folder, IntParam, StrParam, Task, InputSpecs, OutputSpecs,
    TaskInputs, TaskOutputs, task_decorator, InputSpec, OutputSpec, ConfigSpecs
)
from gws_omix import FastqFolder
from .salmon_env import SalmonShellProxyHelper  # A shell proxy for Salmon

@task_decorator("Salmon_Quant", human_name="Salmon_Quant",
                short_description="Quantify RNA-seq reads to a reference transcriptome using Salmon")
class SalmonQuant(Task):
    """
    This task runs Salmon to quantify RNA-seq reads, similar to how
    the Hisat2 script matches filenames. We do:
      - Single-end: any .fastq.gz -> run salmon quant -r
      - Paired-end: we match forward & reverse with user-defined separators.

    For each sample, we create a subfolder <sample_name> in 'result',
    containing the usual Salmon output (quant.sf, logs, etc.).
    """

    input_specs: InputSpecs = InputSpecs({
        'fastq_folder': InputSpec(FastqFolder, human_name="Trimmed Fastq Folder",
                                  short_description="Folder containing trimmed FASTQ reads"),
        'salmon_index': InputSpec(Folder, human_name="Salmon Index",
                                  short_description="Folder containing Salmon transcriptome index")
    })

    output_specs: OutputSpecs = OutputSpecs({
        'output': OutputSpec(Folder, human_name="Salmon Quant Folder",
                             short_description="Folder containing one Salmon quant subfolder per sample")
    })

    config_specs: ConfigSpecs = {
        "threads": IntParam(default_value=4, min_value=2, short_description="Number of threads"),
        "sequencing_type": StrParam(allowed_values=["Paired-end", "Single-end"],
                                    short_description="Sequencing type"),
        "Forward_separator": StrParam(allowed_values=["R1", "1", "r1", " "],
                                      short_description="Forward read identifier for Paired-end"),
        "Reverse_separator": StrParam(allowed_values=["R2", "2", "r2", " "],
                                      short_description="Reverse read identifier for Paired-end")
    }

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Main entry point for Salmon quantification. """
        salmon_index: Folder = inputs['salmon_index']
        fastq_folder: FastqFolder = inputs['fastq_folder']

        threads = params["threads"]
        seq_type = params["sequencing_type"]
        fwd_sep = params["Forward_separator"]
        rev_sep = params["Reverse_separator"]

        # The Salmon index path (e.g. "salmon_index" folder)
        salmon_index_path = salmon_index.path

        # Create output directory
        shell_proxy = SalmonShellProxyHelper.create_proxy(self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'result')
        os.makedirs(result_path, exist_ok=True)

        # Patterns to detect forward / reverse files
        forward_pattern = re.compile(rf"^(.*?)[-_.](?:{re.escape(fwd_sep)})\w*\.fastq\.gz$")
        reverse_pattern = re.compile(rf"^(.*?)[-_.](?:{re.escape(rev_sep)})\w*\.fastq\.gz$")

        if seq_type == "Single-end":
            print("[INFO] Single-end quantification with Salmon...")
            # For single-end, every file .fastq.gz is a separate sample
            for fname in os.listdir(fastq_folder.path):
                if fname.endswith(".fastq.gz"):
                    sample_name = fname.replace(".fastq.gz", "")
                    fastq_path = os.path.join(fastq_folder.path, fname)

                    # Output subfolder
                    out_dir = os.path.join(result_path, sample_name)
                    os.makedirs(out_dir, exist_ok=True)

                    salmon_cmd = (
                        f"salmon quant "
                        f"-i {salmon_index_path} "
                        f"-p {threads} "
                        f"-l A "
                        f"-r {fastq_path} "
                        f"--validateMappings "
                        f"-o {out_dir}"
                    )
                    print("[DEBUG] Single-end CMD:", salmon_cmd)
                    rc = shell_proxy.run(salmon_cmd, shell_mode=True)
                    if rc != 0:
                        raise Exception(f"Error during Salmon quant (single-end) for {sample_name}")

        else:
            print("[INFO] Paired-end quantification with Salmon...")

            # 1) Identify forward & reverse read files
            pairs_dict = {}

            for fname in os.listdir(fastq_folder.path):
                if fname.endswith(".fastq.gz"):
                    fwd_match = forward_pattern.match(fname)
                    if fwd_match:
                        sample_name = fwd_match.group(1)
                        pairs_dict.setdefault(sample_name, {})['FWD'] = fname
                        continue

                    rev_match = reverse_pattern.match(fname)
                    if rev_match:
                        sample_name = rev_match.group(1)
                        pairs_dict.setdefault(sample_name, {})['REV'] = fname
                        continue

            # 2) For each sample in pairs_dict, run Salmon with -1 / -2
            for sample_name, files in pairs_dict.items():
                if 'FWD' in files and 'REV' in files:
                    r1_in = os.path.join(fastq_folder.path, files['FWD'])
                    r2_in = os.path.join(fastq_folder.path, files['REV'])

                    out_dir = os.path.join(result_path, sample_name)
                    os.makedirs(out_dir, exist_ok=True)

                    salmon_cmd = (
                        f"salmon quant "
                        f"-i {salmon_index_path} "
                        f"-p {threads} "
                        f"-l A "
                        f"-1 {r1_in} -2 {r2_in} "
                        f"--validateMappings "
                        f"-o {out_dir}"
                    )
                    print("[DEBUG] Paired-end CMD:", salmon_cmd)

                    rc = shell_proxy.run(salmon_cmd, shell_mode=True)
                    if rc != 0:
                        raise Exception(f"Error during Salmon quant (paired-end) for {sample_name}")

        return {'output': Folder(result_path)}
'''

# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com
'''
import os
import re
import pandas as pd

from gws_core import (
    ConfigParams, Folder, IntParam, StrParam, Task, InputSpecs, OutputSpecs,
    TaskInputs, TaskOutputs, task_decorator, InputSpec, OutputSpec, ConfigSpecs, File
)
from gws_omix import FastqFolder
from .salmon_env import SalmonShellProxyHelper  # A shell proxy for Salmon

@task_decorator("Salmon_Quant", human_name="Salmon_Quant",
                short_description="Quantify RNA-seq reads with Salmon, then merge raw counts into one CSV.")
class SalmonQuant(Task):
    """
    This task:
      1) Detects single-end or paired-end FASTQ data by user config,
      2) Matches forward/reverse reads with regex, removing '_trimmed' from subfolder names,
      3) Runs 'salmon quant' for each sample, storing output in 'result/<sampleName>/quant.sf',
      4) Merges all .quant.sf into a single raw count matrix using 'salmon quantmerge --column numreads',
      5) Post-processes columns to remove '.quant.sf' and leftover '_trimmed',
      6) Writes final 'salmon_merged_counts.csv' in the result folder as CSV.

    Filenames must end with .fastq.gz.
    Single-end => each .fastq.gz is its own sample
    Paired-end => forward/reverse separators via user config (R1/R2, 1/2, etc.).
    """

    input_specs: InputSpecs = InputSpecs({
        'fastq_folder': InputSpec(FastqFolder, human_name="Trimmed FASTQ Folder",
                                  short_description="Folder containing trimmed FASTQ reads"),
        'salmon_index': InputSpec(Folder, human_name="Salmon Index",
                                  short_description="Folder containing Salmon transcriptome index")
    })

    output_specs: OutputSpecs = OutputSpecs({
        'output': OutputSpec(Folder, human_name="Salmon Quant Folder",
                             short_description="Folder containing subfolders for each sample"),
        'raw_counts': OutputSpec(File, human_name="Merged Raw Count CSV",
                                 short_description="Single CSV with gene rows and sample columns from NumReads")
    })

    config_specs: ConfigSpecs = {
        "threads": IntParam(default_value=4, min_value=2,
                            short_description="Number of threads"),
        "sequencing_type": StrParam(allowed_values=["Paired-end", "Single-end"],
                                    short_description="Sequencing type"),
        "Forward_separator": StrParam(allowed_values=["R1", "1", "r1", " "],
                                      short_description="Forward read identifier (for Paired-end)"),
        "Reverse_separator": StrParam(allowed_values=["R2", "2", "r2", " "],
                                      short_description="Reverse read identifier (for Paired-end)")
    }

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        # 1) Retrieve user parameters
        salmon_index: Folder = inputs['salmon_index']
        fastq_folder: FastqFolder = inputs['fastq_folder']
        threads = params["threads"]
        seq_type = params["sequencing_type"]
        fwd_sep = params["Forward_separator"]
        rev_sep = params["Reverse_separator"]

        salmon_index_path = salmon_index.path

        # 2) Create output folder
        shell_proxy = SalmonShellProxyHelper.create_proxy(self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'result')
        os.makedirs(result_path, exist_ok=True)

        # Regex for forward & reverse
        forward_pattern = re.compile(rf"^(.*?)[-_.](?:{re.escape(fwd_sep)})\w*\.fastq\.gz$")
        reverse_pattern = re.compile(rf"^(.*?)[-_.](?:{re.escape(rev_sep)})\w*\.fastq\.gz$")

        # Single-end or Paired-end
        if seq_type == "Single-end":
            print("[INFO] Single-end Salmon quant...")
            for fname in os.listdir(fastq_folder.path):
                if fname.endswith(".fastq.gz"):
                    # Remove .fastq.gz => remove _trimmed
                    sample_name = fname.replace(".fastq.gz", "").replace("_trimmed", "")
                    fastq_path = os.path.join(fastq_folder.path, fname)

                    out_dir = os.path.join(result_path, sample_name)
                    os.makedirs(out_dir, exist_ok=True)

                    salmon_cmd = (
                        f"salmon quant "
                        f"-i {salmon_index_path} "
                        f"-p {threads} "
                        f"-l A "
                        f"-r {fastq_path} "
                        f"--validateMappings "
                        f"-o {out_dir}"
                    )
                    print("[DEBUG] Single-end CMD:", salmon_cmd)
                    rc = shell_proxy.run(salmon_cmd, shell_mode=True)
                    if rc != 0:
                        raise Exception(f"Error during Salmon quant (single-end) for {sample_name}")

        else:
            print("[INFO] Paired-end Salmon quant...")
            pairs_dict = {}
            # Identify forward & reverse
            for fname in os.listdir(fastq_folder.path):
                if fname.endswith(".fastq.gz"):
                    fwd_m = forward_pattern.match(fname)
                    if fwd_m:
                        sample_id = fwd_m.group(1).replace("_trimmed", "")
                        pairs_dict.setdefault(sample_id, {})['FWD'] = fname
                        continue

                    rev_m = reverse_pattern.match(fname)
                    if rev_m:
                        sample_id = rev_m.group(1).replace("_trimmed", "")
                        pairs_dict.setdefault(sample_id, {})['REV'] = fname
                        continue

            # For each sample => run salmon
            for sample_id, files in pairs_dict.items():
                if 'FWD' in files and 'REV' in files:
                    r1_in = os.path.join(fastq_folder.path, files['FWD'])
                    r2_in = os.path.join(fastq_folder.path, files['REV'])

                    out_dir = os.path.join(result_path, sample_id)
                    os.makedirs(out_dir, exist_ok=True)

                    salmon_cmd = (
                        f"salmon quant "
                        f"-i {salmon_index_path} "
                        f"-p {threads} "
                        f"-l A "
                        f"-1 {r1_in} -2 {r2_in} "
                        f"--validateMappings "
                        f"-o {out_dir}"
                    )
                    print("[DEBUG] Paired-end CMD:", salmon_cmd)
                    rc = shell_proxy.run(salmon_cmd, shell_mode=True)
                    if rc != 0:
                        raise Exception(f"Error during Salmon quant (paired-end) for {sample_id}")

        # 3) Gather all quant.sf => build a list for salmon quantmerge
        quant_list = []
        for sample_folder in os.listdir(result_path):
            sample_path = os.path.join(result_path, sample_folder)
            if not os.path.isdir(sample_path):
                continue
            quant_sf = os.path.join(sample_path, "quant.sf")
            if os.path.isfile(quant_sf):
                quant_list.append(quant_sf)

        print("[DEBUG] quant_list:", quant_list)  # <-- print the full list
        if not quant_list:
            raise Exception("No 'quant.sf' files found in 'result' subfolders. Possibly no matching reads or Salmon quant failed.")

        # 4) Run salmon quantmerge --column numreads
        merged_txt = os.path.join(result_path, "temp_merged_numreads.txt")
        merge_cmd = (
            "salmon quantmerge "
            f"--quants {' '.join(quant_list)} "
            f"--column numreads "
            f"-o {merged_txt}"
        )
        print("[DEBUG] Salmon Merge command:", merge_cmd)

        rc_merge = shell_proxy.run(merge_cmd, shell_mode=True)
        if rc_merge != 0:
            raise Exception("Error running 'salmon quantmerge' for raw counts.")

        # 5) Post-process columns => remove .quant.sf or leftover _trimmed => write final CSV
        final_csv = os.path.join(result_path, "salmon_merged_counts.csv")

        if not os.path.isfile(merged_txt):
            raise Exception(f"'salmon quantmerge' did not produce {merged_txt} as expected.")

        df = pd.read_csv(merged_txt, sep='\t')
        old_cols = list(df.columns)

        new_cols = []
        for c in old_cols:
            # e.g. "SRR13978640_trimmed.quant.sf" => "SRR13978640"
            c = c.replace(".quant.sf", "").replace("_trimmed", "")
            new_cols.append(c)
        df.columns = new_cols

        df.to_csv(final_csv, index=False)
        print(f"[INFO] Merged raw count matrix => {final_csv}")

        return {
            'output': Folder(result_path),
            'raw_counts': File(final_csv)
        }
'''

import os
import re
import pandas as pd

from gws_core import (
    ConfigParams, Folder, IntParam, StrParam, Task, InputSpecs, OutputSpecs,
    TaskInputs, TaskOutputs, task_decorator, InputSpec, OutputSpec, ConfigSpecs, File
)
from gws_omix import FastqFolder
from .salmon_env import SalmonShellProxyHelper  # A shell proxy for Salmon

@task_decorator("Salmon_Quant", human_name="Salmon_Quant",
                short_description="Quantify RNA-seq reads with Salmon.")
class SalmonQuant(Task):
    """
    This task:
      1) Detects single-end or paired-end FASTQ data by user config,
      2) Matches forward/reverse reads with regex, removing '_trimmed' from subfolder names,
      3) Runs 'salmon quant' for each sample, storing output in 'result/<sampleName>/quant.sf'.

    Filenames must end with .fastq.gz.
    Single-end => each .fastq.gz is its own sample
    Paired-end => forward/reverse separators via user config (R1/R2, 1/2, etc.).
    This script does NOT perform any merging step.
    """

    input_specs: InputSpecs = InputSpecs({
        'fastq_folder': InputSpec(FastqFolder, human_name="Trimmed FASTQ Folder",
                                  short_description="Folder containing trimmed FASTQ reads"),
        'salmon_index': InputSpec(Folder, human_name="Salmon Index",
                                  short_description="Folder containing Salmon transcriptome index")
    })

    output_specs: OutputSpecs = OutputSpecs({
        'output': OutputSpec(Folder, human_name="Salmon Quant Folder",
                             short_description="Folder containing subfolders for each sample")
        # We remove the 'raw_counts' output for now, since no merging is done.
    })

    config_specs: ConfigSpecs = {
        "threads": IntParam(default_value=4, min_value=2,
                            short_description="Number of threads"),
        "sequencing_type": StrParam(allowed_values=["Paired-end", "Single-end"],
                                    short_description="Sequencing type"),
        "Forward_separator": StrParam(allowed_values=["R1", "1", "r1", " "],
                                      short_description="Forward read identifier (for Paired-end)"),
        "Reverse_separator": StrParam(allowed_values=["R2", "2", "r2", " "],
                                      short_description="Reverse read identifier (for Paired-end)")
    }

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        # 1) Retrieve user parameters
        salmon_index: Folder = inputs['salmon_index']
        fastq_folder: FastqFolder = inputs['fastq_folder']
        threads = params["threads"]
        seq_type = params["sequencing_type"]
        fwd_sep = params["Forward_separator"]
        rev_sep = params["Reverse_separator"]

        salmon_index_path = salmon_index.path

        # 2) Create output folder
        shell_proxy = SalmonShellProxyHelper.create_proxy(self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'result')
        os.makedirs(result_path, exist_ok=True)

        # Regex for forward & reverse
        forward_pattern = re.compile(rf"^(.*?)[-_.](?:{re.escape(fwd_sep)})\w*\.fastq\.gz$")
        reverse_pattern = re.compile(rf"^(.*?)[-_.](?:{re.escape(rev_sep)})\w*\.fastq\.gz$")

        # Single-end or Paired-end
        if seq_type == "Single-end":
            print("[INFO] Single-end Salmon quant...")
            for fname in os.listdir(fastq_folder.path):
                if fname.endswith(".fastq.gz"):
                    # Remove .fastq.gz => remove _trimmed
                    sample_name = fname.replace(".fastq.gz", "").replace("_trimmed", "")
                    fastq_path = os.path.join(fastq_folder.path, fname)

                    out_dir = os.path.join(result_path, sample_name)
                    os.makedirs(out_dir, exist_ok=True)

                    salmon_cmd = (
                        f"salmon quant "
                        f"-i {salmon_index_path} "
                        f"-p {threads} "
                        f"-l A "
                        f"-r {fastq_path} "
                        f"--validateMappings "
                        f"-o {out_dir}"
                    )
                    print("[DEBUG] Single-end CMD:", salmon_cmd)
                    rc = shell_proxy.run(salmon_cmd, shell_mode=True)
                    if rc != 0:
                        raise Exception(f"Error during Salmon quant (single-end) for {sample_name}")

        else:
            print("[INFO] Paired-end Salmon quant...")
            pairs_dict = {}

            # Identify forward & reverse
            for fname in os.listdir(fastq_folder.path):
                if fname.endswith(".fastq.gz"):
                    fwd_m = forward_pattern.match(fname)
                    if fwd_m:
                        sample_id = fwd_m.group(1).replace("_trimmed", "")
                        pairs_dict.setdefault(sample_id, {})['FWD'] = fname
                        continue

                    rev_m = reverse_pattern.match(fname)
                    if rev_m:
                        sample_id = rev_m.group(1).replace("_trimmed", "")
                        pairs_dict.setdefault(sample_id, {})['REV'] = fname
                        continue

            # For each sample => run salmon
            for sample_id, files in pairs_dict.items():
                if 'FWD' in files and 'REV' in files:
                    r1_in = os.path.join(fastq_folder.path, files['FWD'])
                    r2_in = os.path.join(fastq_folder.path, files['REV'])

                    out_dir = os.path.join(result_path, sample_id)
                    os.makedirs(out_dir, exist_ok=True)

                    salmon_cmd = (
                        f"salmon quant "
                        f"-i {salmon_index_path} "
                        f"-p {threads} "
                        f"-l A "
                        f"-1 {r1_in} -2 {r2_in} "
                        f"--validateMappings "
                        f"-o {out_dir}"
                    )
                    print("[DEBUG] Paired-end CMD:", salmon_cmd)
                    rc = shell_proxy.run(salmon_cmd, shell_mode=True)
                    if rc != 0:
                        raise Exception(f"Error during Salmon quant (paired-end) for {sample_id}")

        # Return only the folder. No merges are performed.
        return {
            'output': Folder(result_path)
        }


