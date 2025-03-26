import os
import re
from gws_core import (ConfigParams, Folder, InputSpec, OutputSpec, InputSpecs, OutputSpecs, ConfigSpecs,
                      Task, TaskInputs, TaskOutputs, task_decorator, ConfigSpecs, StrParam)

from .multiqc_env import MultiQcShellProxyHelper


@task_decorator("MultiQC", human_name="MultiQC",
                short_description="Create aggregate bioinformatics analysis reports across many samples")
class MultiQc(Task):
    """
    MultiQC is a tool that allows the aggregation and summarization of results from multiple bioinformatics analyses in a single report.
    It simplifies the task of analyzing large datasets by automatically generating interactive and visually appealing reports.
    """

    input_specs = InputSpecs({
        'fastqc_reports_folder': InputSpec(
            Folder, human_name="Fastqc Reports Files", short_description="Output statistics files of QC tool")
    })

    output_specs = OutputSpecs({
        'output': OutputSpec(Folder, human_name="Combined Quality Report",
                             short_description="Quality report for all the fastq files")
    })

    config_specs = ConfigSpecs({
        "sequencing_mode": StrParam(default_value="paired", allowed_values=["paired", "single"],
                                    short_description="Sequencing mode: paired-end or single-end")
    })

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # Récupérer les paramètres d'entrée
        input_folder: Folder = inputs['fastqc_reports_folder']
        sequencing_mode = params["sequencing_mode"]

        # Créer le répertoire de travail
        shell_proxy = MultiQcShellProxyHelper.create_proxy(
            self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'result')
        os.makedirs(result_path, exist_ok=True)

        # Définition du pattern pour identifier les fichiers FastQC (html et zip)
        qc_pattern = re.compile(r'.*\.(html|zip)$')

        # Récupérer tous les fichiers FastQC disponibles
        fastqc_files = [f for f in os.listdir(
            input_folder.path) if qc_pattern.match(f)]

        if sequencing_mode == "single":
            # Vérifier qu'on a bien trouvé des fichiers
            if not fastqc_files:
                raise Exception(
                    "Aucun rapport FastQC trouvé pour les fichiers single-end.")

            fastqc_cmd = f"multiqc {' '.join([os.path.join(input_folder.path, f) for f in fastqc_files])} -o {result_path} -n multiqc_report.html"
            res = shell_proxy.run(fastqc_cmd, shell_mode=True)

        else:  # Mode paired-end
            # Patterns pour R1 et R2
            r1_pattern = re.compile(r'.*?[-_.][R1|1|r1]\w*\.(html|zip)')
            r2_pattern = re.compile(r'.*?[-_.][R2|2|r2]\w*\.(html|zip)')

            r1_files = [f for f in fastqc_files if r1_pattern.match(f)]
            r2_files = [f for f in fastqc_files if r2_pattern.match(f)]

            # Vérification de la présence des fichiers R1 et R2
            if not r1_files or not r2_files:
                raise Exception(
                    "Aucun rapport FastQC trouvé pour les fichiers paired-end.")

            forward_path = os.path.join(result_path, 'forward')
            reverse_path = os.path.join(result_path, 'reverse')
            os.makedirs(forward_path, exist_ok=True)
            os.makedirs(reverse_path, exist_ok=True)

            # Commandes MultiQC
            command_forward = f"multiqc {' '.join([os.path.join(input_folder.path, f) for f in r1_files])} -o {forward_path} -n multiqc_forward.html"
            command_reverse = f"multiqc {' '.join([os.path.join(input_folder.path, f) for f in r2_files])} -o {reverse_path} -n multiqc_reverse.html"
            command_combined = f"multiqc {' '.join([os.path.join(input_folder.path, f) for f in fastqc_files])} -o {result_path} -n multiqc_combined.html"

            # Exécution des commandes
            res_forward = shell_proxy.run(command_forward, shell_mode=True)
            res_reverse = shell_proxy.run(command_reverse, shell_mode=True)
            res_combined = shell_proxy.run(command_combined, shell_mode=True)

            if res_forward != 0 or res_reverse != 0 or res_combined != 0:
                raise Exception(
                    "Une erreur est survenue lors de la génération des rapports MultiQC.")

        # Retourner le dossier de résultats
        return {'output': Folder(result_path)}
