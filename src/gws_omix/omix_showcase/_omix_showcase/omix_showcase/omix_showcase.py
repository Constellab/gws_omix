from typing import Optional

import pandas as pd
import reflex as rx
from gws_core import File
from gws_reflex_base import get_theme, render_main_container
from gws_reflex_main import ReflexMainState


class State(ReflexMainState):

    @rx.var
    def database_data(self) -> pd.DataFrame:
        resources = self.get_resources()
        if not resources or len(resources) == 0:
            return pd.DataFrame()
        data_file: File = resources[0]
        if not data_file:
            return pd.DataFrame()
        data: pd.DataFrame = pd.read_excel(data_file.path, sheet_name=0)
        if data is None or data.empty:
            return pd.DataFrame()
        return data

    @rx.var
    def constellab_data(self) -> pd.DataFrame:
        resources = self.get_resources()
        if not resources or len(resources) == 0:
            return pd.DataFrame()
        data_file: File = resources[0]
        if not data_file:
            return pd.DataFrame()
        data: pd.DataFrame = pd.read_excel(data_file.path, sheet_name=1)
        if data is None or data.empty:
            return pd.DataFrame()
        return data


app = rx.App(
    theme=get_theme(),
    stylesheets=[
        "/style.css"
    ]
)


@rx.page(title="Functional annotation through orthology assignment standalone application", on_load=ReflexMainState.on_load)
def index():

    # Render the main container with the app content.
    # The content will be displayed once the state is initialized.
    # If the state is not initialized, a loading spinner will be shown.
    return render_main_container(
        rx.flex(
            rx.heading("Functional annotation through orthology assignment standalone application", font_size="2em"),
            rx.tabs.root(
                rx.tabs.list(
                    rx.tabs.trigger("Introduction", value="introduction", class_name="clickable"),
                    rx.tabs.trigger("Database Table", value="database_table", class_name="clickable"),
                    rx.tabs.trigger("Constellab Table", value="constellab_table", class_name="clickable"),),
                rx.tabs.content(
                    rx.box(
                        rx.text(
                            "EggNOG-mapper is a powerful tool for the fast and precise functional annotation of novel sequences such as proteins, coding sequences (CDS), genomes, and metagenomes. It leverages precomputed orthologous groups (OGs) and phylogenies from the ",
                            rx.link(
                                "EggNOG Database", href="http://eggnog5.embl.de/download/emapperdb-5.0.2/",
                                target="_blank"),
                            " to ensure annotations are transferred only from fine grained orthologs, avoiding misleading annotations from close paralogs.",),
                        rx.text(
                            "This tool significantly improves annotation accuracy compared to traditional homology-based tools like BLAST, making it a preferred solution for researchers working on novel genomes, transcriptomes, or metagenomic datasets."),
                        rx.text(
                            "EggNOG-mapper has been wrapped into a Constellab Task, enabling easy integration into bioinformatics workflows with reproducibility and scalability."),
                        rx.text(
                            "The input file MGYG000307600.faa was downloaded from the ", rx.link(
                                "mgnify database",
                                href="https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/chicken-gut/v1.0.1/species_catalogue/MGYG0003076/MGYG000307600/genome/",
                                target="_blank"),
                            "."),
                        rx.text(
                            "This FASTA file was used as input for the EggNOG-mapper pipeline. After running the script, the resulting annotation file was compared to the functional annotations published in the mgnify database."),
                        rx.text(
                            "The results were consistent and comparable, demonstrating the reliability and accuracy of the pipeline for functional annotation."),
                        rx.text(
                            "The documentation of the task is available in Constellab community ", rx.link(
                                "here",
                                href="https://constellab.community/bricks/gws_omix/latest/doc/use-cases/functional-annotation-through-orthology-assignment/ac867645-9bcf-4f36-9893-b3da27431e16",
                                target="_blank"),
                            "."),
                        rx.button(
                            "Download FAA Input File", on_click=rx.download(
                                url=rx.asset("MGYG000307600.faa.txt"),
                                filename="MGYG000307600.faa.txt"),
                            id="download_faa_file_button", style={"margin-top": "1em"}, class_name="clickable"),
                        class_name="introduction-box"),
                    value="introduction", class_name="tab-content"),
                rx.tabs.content(
                    rx.data_table(
                        data=State.database_data, pagination=True, resizable=True, width="100%",
                        class_name="data-table"),
                    value="database_table", class_name="tab-content"),
                rx.tabs.content(
                    rx.data_table(
                        data=State.constellab_data, pagination=True, resizable=True, width="100%",
                        class_name="data-table"),
                    value="constellab_table", class_name="tab-content"),
                default_value="introduction"),
            direction="column", width="100%", spacing="2", padding="1em"))
