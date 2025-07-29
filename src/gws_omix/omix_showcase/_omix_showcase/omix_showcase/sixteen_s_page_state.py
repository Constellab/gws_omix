import json
from typing import Any, Dict, Optional

import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import reflex as rx
from gws_core import File, Folder
from gws_reflex_main import ReflexMainState


class SixteenSPageState(ReflexMainState):

    @rx.var(cache=True)
    def quality_check_forward_reads_url(self) -> str:
        """Return the URL for the quality check forward reads."""
        return self.get_param("quality_check_forward_reads", "")

    @rx.var(cache=True)
    def quality_check_reverse_reads_url(self) -> str:
        """Return the URL for the quality check reverse reads."""
        return self.get_param("quality_check_reverse_reads", "")

    @rx.var(cache=True)
    def feature_inference(self) -> pd.DataFrame:
        sixteen_s_resource_folder = self._get_sixteen_s_resource_folder()
        if not sixteen_s_resource_folder:
            return pd.DataFrame()
        feature_inference_file: File = sixteen_s_resource_folder.get_sub_node("feature_inference.csv")
        return self._check_and_return_csv_file(feature_inference_file)

    @rx.var(cache=True)
    def rarefaction_analysis_observed_features_url(self) -> Optional[str]:
        """Return the URL for the rarefaction analysis observed features."""
        return self.get_param("rarefaction_analysis_observed_features")

    @rx.var(cache=True)
    def rarefaction_analysis_shannon_index_url(self) -> Optional[str]:
        """Return the URL for the rarefaction analysis shannon index."""
        return self.get_param("rarefaction_analysis_shannon_index")

    @rx.var(cache=True)
    def diversity_analysis(self) -> dict[str, dict[str, pd.DataFrame]]:
        res = {}
        sixteen_s_resource_folder = self._get_sixteen_s_resource_folder()
        if not sixteen_s_resource_folder:
            return res
        diversity_analysis_folder: Folder = sixteen_s_resource_folder.get_sub_node('diversity_analysis')
        if not diversity_analysis_folder:
            return res

        file_names = diversity_analysis_folder.list_dir()
        if not file_names:
            return rx.text("No files found in the diversity analysis folder.")

        alpha_diversity_files = [
            file for file in file_names if file.startswith("Alpha") and file.endswith(".csv")
        ]

        beta_diversity_files = [
            file for file in file_names if file.startswith("Beta") and file.endswith(".csv")
        ]

        if not alpha_diversity_files and not beta_diversity_files:
            return res

        if alpha_diversity_files:
            res['Alpha Diversity'] = {}
            for file_name in alpha_diversity_files:
                file: File = diversity_analysis_folder.get_sub_node(file_name)
                if file and isinstance(file, File):
                    res['Alpha Diversity'][file_name] = self._check_and_return_csv_file(file)

        if beta_diversity_files:
            res['Beta Diversity'] = {}
            for file_name in beta_diversity_files:
                file: File = diversity_analysis_folder.get_sub_node(file_name)
                if file and isinstance(file, File):
                    res['Beta Diversity'][file_name] = self._check_and_return_csv_file(file)
        return res

    @rx.var(cache=True)
    def get_taxonomy_analysis_kingdom_url(self) -> Optional[str]:
        """Return the URL for the taxonomy analysis kingdom."""
        return self.get_param("taxonomy_analysis_kingdom")

    @rx.var(cache=True)
    def get_taxonomy_analysis_phylum_url(self) -> Optional[str]:
        """Return the URL for the taxonomy analysis phylum."""
        return self.get_param("taxonomy_analysis_phylum")

    @rx.var(cache=True)
    def get_taxonomy_analysis_class_url(self) -> Optional[str]:
        """Return the URL for the taxonomy analysis class."""
        return self.get_param("taxonomy_analysis_class")

    @rx.var(cache=True)
    def get_taxonomy_analysis_order_url(self) -> Optional[str]:
        """Return the URL for the taxonomy analysis order."""
        return self.get_param("taxonomy_analysis_order")

    @rx.var(cache=True)
    def get_taxonomy_analysis_family_url(self) -> Optional[str]:
        """Return the URL for the taxonomy analysis family."""
        return self.get_param("taxonomy_analysis_family")

    @rx.var(cache=True)
    def get_taxonomy_analysis_genus_url(self) -> Optional[str]:
        """Return the URL for the taxonomy analysis genus."""
        return self.get_param("taxonomy_analysis_genus")

    @rx.var(cache=True)
    def get_taxonomy_analysis_species_url(self) -> Optional[str]:
        """Return the URL for the taxonomy analysis species."""
        return self.get_param("taxonomy_analysis_species")

    @rx.var(cache=True)
    def get_functional_analysis_plot(self) -> Optional[go.Figure]:
        sixteen_s_resource_folder = self._get_sixteen_s_resource_folder()
        if not sixteen_s_resource_folder:
            return None
        functional_analysis_folder: Folder = sixteen_s_resource_folder.get_sub_node("functional_analysis_prediction")
        if not functional_analysis_folder:
            return None
        functional_analysis_file_plot: File = functional_analysis_folder.get_sub_node("plotly.json")
        if not functional_analysis_file_plot or not isinstance(functional_analysis_file_plot, File):
            return None
        fig = pio.from_json(functional_analysis_file_plot.read())
        return fig

    @rx.var(cache=True)
    def get_functional_analysis_data(self) -> pd.DataFrame:
        functional_analysis_folder_resource_set: Folder = self._get_functional_analysis_folder_resource_set()
        if not functional_analysis_folder_resource_set:
            return pd.DataFrame()
        functional_analysis_file_data: File = functional_analysis_folder_resource_set.get_sub_node(
            "daa_annotated_results_HC.csv")
        if not functional_analysis_file_data or not isinstance(functional_analysis_file_data, File):
            return pd.DataFrame()
        data: pd.DataFrame = pd.read_csv(functional_analysis_file_data.path)
        if data is None or data.empty:
            return pd.DataFrame()
        return data

    @rx.var(cache=True)
    def get_functional_analysis_pathway_errorbar_image(self) -> Optional[str]:
        functional_analysis_folder_resource_set: Folder = self._get_functional_analysis_folder_resource_set()
        if not functional_analysis_folder_resource_set:
            return None
        functional_analysis_file_pathway_errorbar_image: File = functional_analysis_folder_resource_set.get_sub_node(
            "pathway_errorbar_HC.png")
        if not functional_analysis_file_pathway_errorbar_image or not isinstance(
                functional_analysis_file_pathway_errorbar_image, File):
            return None
        return functional_analysis_file_pathway_errorbar_image.path

    @rx.var(cache=True)
    def get_functional_analysis_pathway_heatmap_image(self) -> Optional[str]:
        functional_analysis_folder_resource_set: Folder = self._get_functional_analysis_folder_resource_set()
        if not functional_analysis_folder_resource_set:
            return None
        functional_analysis_file_pathway_heatmap_image: File = functional_analysis_folder_resource_set.get_sub_node(
            "pathway_errorbar_HC.png")
        if not functional_analysis_file_pathway_heatmap_image or not isinstance(
                functional_analysis_file_pathway_heatmap_image, File):
            return None
        return functional_analysis_file_pathway_heatmap_image.path

    def _get_sixteen_s_resource_folder(self) -> Optional[Folder]:
        resources = self.get_resources()
        if not resources or len(resources) != 2:
            return None
        sixteen_s_resource_folder: Folder = resources[1]
        if not sixteen_s_resource_folder:
            return None
        return sixteen_s_resource_folder

    def _get_quality_check_data_file(self, file_name: str) -> Optional[File]:
        sixteen_s_resource_folder = self._get_sixteen_s_resource_folder()
        if not sixteen_s_resource_folder:
            return None
        quality_check_folder: Folder = sixteen_s_resource_folder.get_sub_node('quality_check')
        if not quality_check_folder:
            return None
        data_file: File = quality_check_folder.get_sub_node(f"{file_name}")
        if not data_file or not isinstance(data_file, File):
            return None
        return data_file

    def _get_rarefaction_analysis_file(self, csv_file_name: str) -> Optional[File]:
        sixteen_s_resource_folder = self._get_sixteen_s_resource_folder()
        if not sixteen_s_resource_folder:
            return None
        rarefaction_analysis_folder: Folder = sixteen_s_resource_folder.get_sub_node('rarefaction_analysis')
        if not rarefaction_analysis_folder:
            return None
        data_file: File = rarefaction_analysis_folder.get_sub_node(f"{csv_file_name}.csv")
        if not data_file or not isinstance(data_file, File):
            return None
        return data_file

    def _check_and_return_csv_file(self, file: File) -> pd.DataFrame:
        if not file or not isinstance(file, File):
            return pd.DataFrame()
        data: pd.DataFrame = pd.read_csv(file.path)
        if data is None or data.empty:
            return pd.DataFrame()
        return data

    def _get_functional_analysis_folder_resource_set(self) -> Optional[Folder]:
        sixteen_s_resource_folder = self._get_sixteen_s_resource_folder()
        if not sixteen_s_resource_folder:
            return None
        functional_analysis_folder: Folder = sixteen_s_resource_folder.get_sub_node("functional_analysis_prediction")
        if not functional_analysis_folder:
            return None
        functional_analysis_folder_resource_set: Folder = functional_analysis_folder.get_sub_node("resource_set")
        if not functional_analysis_folder_resource_set:
            return None
        return functional_analysis_folder_resource_set
