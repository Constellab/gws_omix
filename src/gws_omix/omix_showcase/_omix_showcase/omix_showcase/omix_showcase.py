
from typing import List

import reflex as rx
from gws_omix.omix_showcase._omix_showcase.omix_showcase.eggnog_page_state import \
    EggnogPageState
from gws_omix.omix_showcase._omix_showcase.omix_showcase.sixteen_s_page_state import \
    SixteenSPageState
from gws_reflex_base import get_theme, render_main_container
from gws_reflex_main import ReflexMainState


def download_link_button(label: str, href: str) -> rx.Component:
    return rx.link(
        label,
        href=href,
        class_name="download-href-button"
    )


def sidebar_item(
    text: str, icon: str, href: str
) -> rx.Component:
    return rx.link(
        rx.hstack(
            rx.icon(icon, class_name="sidebar-icon"),
            rx.text(text, size="3"),
            width="100%",
            padding_x="0.5rem",
            padding_y="0.75rem",
            align="center",
            style={
                "_hover": {
                    "bg": rx.color("accent", 4),
                    "color": rx.color("accent", 11),
                },
                "border-radius": "0.5em",
            },
        ),
        href=href,
        underline="none",
        weight="medium",
        width="100%",
    )


def sidebar_items() -> rx.Component:
    return rx.vstack(
        sidebar_item("Home", "home", "/#"),
        sidebar_item(
            "Eggnog", "square-library",
            "/eggnog"),
        sidebar_item("16S RNA", "square-library", "/16s-rna"),
        spacing="1", width="100%",)


def sidebar() -> rx.Component:
    return rx.box(
        rx.desktop_only(
            rx.vstack(
                rx.hstack(
                    rx.image(
                        src=rx.asset("omix_logo.png"),
                        width="2.25em",
                        height="auto",
                        border_radius="25%",
                    ),
                    rx.heading(
                        "Omix Showcases", size="7", weight="bold"
                    ),
                    align="center",
                    justify="start",
                    padding_x="0.5rem",
                    width="100%",
                ),
                sidebar_items(),
                spacing="5",
                # position="fixed",
                # left="0px",
                # top="0px",
                # z_index="5",
                padding_x="1em",
                padding_y="1.5em",
                bg=rx.color("accent", 3),
                align="start",
                # height="100%",
                height="650px",
                width="16em",
                class_name="sidebar-container"),
            class_name="sidebar-container"),
        rx.mobile_and_tablet(
            rx.drawer.root(
                rx.drawer.trigger(
                    rx.icon("align-justify", size=30)
                ),
                rx.drawer.overlay(z_index="5"),
                rx.drawer.portal(
                    rx.drawer.content(
                        rx.vstack(
                            rx.box(
                                rx.drawer.close(
                                    rx.icon("x", size=30)
                                ),
                                width="100%",
                            ),
                            sidebar_items(),
                            spacing="5",
                            width="100%",
                        ),
                        top="auto",
                        right="auto",
                        height="100%",
                        width="20em",
                        padding="1.5em",
                        bg=rx.color("accent", 2),
                    ),
                    width="100%",
                ),
                direction="left",
            ),
            padding="1em",
        ),
        class_name="sidebar-container")


def page_component(page_box_content: rx.Component) -> rx.Component:
    return rx.flex(
        sidebar(),
        rx.box(page_box_content, class_name="page_container")
    )


app = rx.App(
    theme=get_theme(),
    stylesheets=[
        "/style.css"
    ]


)


@rx.page(title="Home", on_load=ReflexMainState.on_load)
def index():
    # Render the main container with the app content.
    # The content will be displayed once the state is initialized.
    # If the state is not initialized, a loading spinner will be shown.
    return render_main_container(
        page_component(
            rx.flex(
                rx.heading("Omix Showcase", font_size="2em"),
                rx.text("Welcome to the Omix Showcase!"),
                rx.text("This application demonstrates various functionalities of the Omix platform."),
                class_name="overview-box",
                width="100%",
                direction="column",
                spacing="2",
                padding="1em")))


@rx.page(title="Eggnog",
         on_load=ReflexMainState.on_load, route="/eggnog")
def functional_annotation_through_orthology_assignment_page():

    # Render the main container with the app content.
    # The content will be displayed once the state is initialized.
    # If the state is not initialized, a loading spinner will be shown.
    return render_main_container(
        page_component(
            rx.flex(
                rx.heading(
                    "Functional annotation through orthology assignment standalone application", font_size="2em"),
                rx.tabs.root(
                    rx.tabs.list(
                        rx.tabs.trigger("Overview", value="overview", class_name=["clickable", "tab-button"]),
                        rx.tabs.trigger("Database Table", value="database_table", class_name=["clickable", "tab-button"]),
                        rx.tabs.trigger("Constellab Table", value="constellab_table", class_name=["clickable", "tab-button"]),),
                    rx.tabs.content(
                        rx.box(
                            rx.heading("Introduction", font_size="1.3em"),
                            rx.text(
                                rx.text.strong("EggNOG-mapper"),
                                " is a powerful tool for the ",
                                rx.text.strong("fast and precise functional annotation"),
                                " of novel sequences such as proteins, coding sequences (CDS), genomes, and metagenomes. It leverages ",
                                rx.text.strong("precomputed orthologous groups (OGs)"),
                                " and ", rx.text.strong("phylogenies"),
                                " from the ", rx.link(
                                    "EggNOG Database", href="http://eggnog5.embl.de/download/emapperdb-5.0.2/",
                                    target="_blank"),
                                " to ensure annotations are transferred only from ", rx.text.strong(
                                    "fine grained orthologs"),
                                ", avoiding misleading annotations from close paralogs.",),
                            rx.text(
                                "This tool significantly improves annotation accuracy compared to traditional homology-based tools like BLAST, making it a preferred solution for researchers working on novel genomes, transcriptomes, or metagenomic datasets."),
                            rx.text(
                                "EggNOG-mapper has been wrapped into a Constellab Task, enabling easy integration into bioinformatics workflows with reproducibility and scalability."),
                            rx.heading("Input file", font_size="1.3em"),
                            rx.flex(
                                rx.text(
                                    "The input file MGYG000307600.faa was downloaded from the ", rx.link(
                                        "mgnify database",
                                        href="https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/chicken-gut/v1.0.1/species_catalogue/MGYG0003076/MGYG000307600/genome/",
                                        target="_blank"),
                                    "."),
                                rx.button(
                                    "Download FAA Input File", on_click=rx.download(
                                        url=rx.asset("MGYG000307600.faa.txt"),
                                        filename="MGYG000307600.faa.txt"),
                                    id="download_faa_file_button", class_name="clickable"),
                                spacing="2", align="center"),
                            rx.heading("Validation with Publicly Available Dataset", font_size="1.3em"),
                            rx.text(
                                "This FASTA file was used as input for the EggNOG-mapper pipeline. After running the script, the resulting annotation file was compared to the functional annotations published in the mgnify database."),
                            rx.text(
                                "The results were consistent and comparable, demonstrating the reliability and accuracy of the pipeline for functional annotation."),
                            rx.heading("Documentation", font_size="1.3em"),
                            rx.text(
                                "The documentation of the task is available in Constellab community ", rx.link(
                                    "here",
                                    href="https://constellab.community/bricks/gws_omix/latest/doc/use-cases/functional-annotation-through-orthology-assignment/ac867645-9bcf-4f36-9893-b3da27431e16",
                                    target="_blank"),
                                "."),
                            class_name="overview-box"),
                        value="overview", class_name="tab-content"),
                    rx.tabs.content(
                        rx.data_table(
                            data=EggnogPageState.database_data, pagination=True, resizable=True, width="100%",
                            class_name="data-table"),
                        download_link_button(
                            "Download Database Table", rx.get_upload_url(EggnogPageState.download_database_data)),
                        value="database_table", class_name="tab-content"),
                    rx.tabs.content(
                        rx.data_table(
                            data=EggnogPageState.constellab_data, pagination=True, resizable=True, width="100%",
                            class_name="data-table"),
                        download_link_button(
                            "Download Constellab Table", rx.get_upload_url(
                                EggnogPageState.download_constellab_data)),
                        value="constellab_table", class_name="tab-content"),
                    default_value="overview"),
                direction="column", width="100%", spacing="2", padding="1em")))


def render_diversity_df(df: List) -> rx.Component:
    return rx.vstack(
        rx.text(df[0]),
        rx.data_table(
            data=df[1],
            pagination=True,
            width="100%", class_name="data-table"))


def diversity_analysis_component(diversity_section: List) -> rx.Component:
    return rx.vstack(
        rx.heading(diversity_section[0], font_size="1.2em"),
        rx.foreach(
            diversity_section[1],
            render_diversity_df
        ),
        class_name="diversity-section"
    )


@rx.page(title="16s RNA", on_load=ReflexMainState.on_load, route="/16s-rna")
def sixteen_s_page():
    # Render the main container with the app content.
    # The content will be displayed once the state is initialized.
    # If the state is not initialized, a loading spinner will be shown.

    return render_main_container(
        page_component(
            rx.flex(
                rx.heading("16S rRNA short-read sequencing analysis", font_size="2em"),
                rx.tabs.root(
                    rx.tabs.list(
                        rx.tabs.trigger("Overview", value="overview", class_name=["clickable", "tab-button"]),
                        rx.tabs.trigger("Quality check", value="quality-check", class_name=["clickable", "tab-button"]),
                        rx.tabs.trigger(
                            "Feature inference", value="feature-inference", class_name=["clickable", "tab-button"]),
                        rx.tabs.trigger(
                            "Rarefaction analysis", value="rarefaction-analysis", class_name=["clickable", "tab-button"]),
                        rx.tabs.trigger(
                            "Diversity analysis", value="diversity-analysis", class_name=["clickable", "tab-button"]),
                        rx.tabs.trigger(
                            "Taxonomy analysis", value="taxonomy-analysis", class_name=["clickable", "tab-button"]),
                        rx.tabs.trigger(
                            "Functional analysis prediction", value="functional-analysis-prediction",
                            class_name=["clickable", "tab-button"]),),
                    rx.tabs.content(
                        rx.box(
                            rx.heading("Introduction", font_size="1.3em"),
                            rx.text(
                                "In 2019, ", rx.text.em("Pietrucci et al."),
                                " published a study on gut microbiota dysbiosis in Parkinson's patients, analyzing 152 fecal samples from 80 patients and 72 healthy controls. The authors conducted 16S ribosomal RNA gene amplicon sequencing combined with dietary/lifestyle data. Sequencing files are available in Sequence Read Archive (PRJNA510730). Their analysis found significantly higher levels of Lactobacillaceae, Enterobacteriaceae, and Enterococcaceae families, and reduced levels of Lachnospiraceae in PD patients compared to controls. Constellab's gws_ubiome brick offers metabarcoding analysis using Qiime2 and Picrust2. Users can upload fastq files and metadata, deploy the pipeline, adjust parameters, and access results in minutes.",
                                as_="div"),
                            rx.heading("Application", font_size="1.3em"),
                            rx.text(
                                "This application is designed exclusively for data visualization. Within it, you can explore read quality metrics, view a global overview of reads after denoising, examine taxonomy classification from kingdom to species, observe rarefaction curves, and explore functional prediction using PICRUSt2. Comprehensive documentation for this module is available on ",
                                rx.link(
                                    "Constellab Community", href=SixteenSPageState.quality_check_forward_reads_url,
                                    target="_blank"),
                                ". We also encourage you to explore the ", rx.link(
                                    "FAQ section", href=SixteenSPageState.quality_check_reverse_reads_url,
                                    target="_blank"),
                                as_="div"),
                            rx.heading("Use Case", font_size="1.3em"),
                            rx.text(
                                "A use case demonstrating the steps and tasks involved in processing the data from this article is available in Constellab community ",
                                rx.link(
                                    "here",
                                    href="https://constellab.community/stories/bf3da2c6-4475-4289-aa7c-b486b23a8bc4/16s-data-analysis-pipeline",
                                    target="_blank"),
                                ".", as_="div"),
                            class_name="overview-box"),
                        value="overview", class_name="tab-content"),
                    rx.tabs.content(
                        rx.heading("Quality Check", font_size="1.3em"),
                        rx.heading("Forward Reads", font_size="1.2em"),
                        rx.el.iframe(
                            src=SixteenSPageState.quality_check_forward_reads_url,
                            width="100%", height="600px", class_name="quality-check-iframe"),
                        rx.heading("Reverse Reads", font_size="1.2em"),
                        rx.el.iframe(
                            src=SixteenSPageState.quality_check_reverse_reads_url,
                            width="100%", height="600px", class_name="quality-check-iframe"),
                        value="quality-check", class_name="tab-content"),
                    rx.tabs.content(
                        rx.heading("Feature Inference", font_size="1.3em"),
                        rx.data_table(
                            data=SixteenSPageState.feature_inference, pagination=True, resizable=True, width="100%",
                            class_name="data-table"),
                        value="feature-inference", class_name="tab-content"),
                    rx.tabs.content(
                        rx.heading("Rarefaction Analysis", font_size="1.3em"),
                        rx.heading("Observed features", font_size="1.2em"),
                        rx.el.iframe(
                            src=SixteenSPageState.rarefaction_analysis_observed_features_url,
                            width="100%", height="600px", class_name="quality-check-iframe"),
                        rx.heading("Shannon index", font_size="1.2em"),
                        rx.el.iframe(
                            src=SixteenSPageState.rarefaction_analysis_shannon_index_url,
                            width="100%", height="600px", class_name="quality-check-iframe"),
                        value="rarefaction-analysis", class_name="tab-content"),
                    rx.tabs.content(
                        rx.heading("Diversity Analysis", font_size="1.3em"),
                        rx.foreach(SixteenSPageState.diversity_analysis, diversity_analysis_component),
                        value="diversity-analysis", class_name="tab-content"),
                    rx.tabs.content(
                        rx.heading("Taxonomy Analysis", font_size="1.3em"),
                        rx.heading("Kingdom", font_size="1.2em"),
                        rx.el.iframe(
                            src=SixteenSPageState.get_taxonomy_analysis_kingdom_url,
                            width="100%", height="600px", class_name="quality-check-iframe"),
                        rx.heading("Phylum", font_size="1.2em"),
                        rx.el.iframe(
                            src=SixteenSPageState.get_taxonomy_analysis_phylum_url,
                            width="100%", height="600px", class_name="quality-check-iframe"),
                        rx.heading("Class", font_size="1.2em"),
                        rx.el.iframe(
                            src=SixteenSPageState.get_taxonomy_analysis_class_url,
                            width="100%", height="600px", class_name="quality-check-iframe"),
                        rx.heading("Order", font_size="1.2em"),
                        rx.el.iframe(
                            src=SixteenSPageState.get_taxonomy_analysis_order_url,
                            width="100%", height="600px", class_name="quality-check-iframe"),
                        rx.heading("Family", font_size="1.2em"),
                        rx.el.iframe(
                            src=SixteenSPageState.get_taxonomy_analysis_family_url,
                            width="100%", height="600px", class_name="quality-check-iframe"),
                        rx.heading("Genus", font_size="1.2em"),
                        rx.el.iframe(
                            src=SixteenSPageState.get_taxonomy_analysis_genus_url,
                            width="100%", height="600px", class_name="quality-check-iframe"),
                        rx.heading("Species", font_size="1.2em"),
                        rx.el.iframe(
                            src=SixteenSPageState.get_taxonomy_analysis_species_url,
                            width="100%", height="600px", class_name="quality-check-iframe"),
                        value="taxonomy-analysis", class_name="tab-content"),
                    rx.tabs.content(
                        rx.heading("Functional Analysis Prediction", font_size="1.3em"),
                        rx.heading("Functional Analysis Plot", font_size="1.2em"),
                        rx.plotly(data=SixteenSPageState.get_functional_analysis_plot, width="100%"),
                        rx.heading("Functional Analysis Data", font_size="1.2em"),
                        rx.data_table(
                            data=SixteenSPageState.get_functional_analysis_data, pagination=True, resizable=True,
                            width="100%", class_name="data-table"),
                        rx.heading("Pathway Heatmap", font_size="1.2em"),
                        rx.image(
                            src=rx.asset("functional_analysis_pathway/pathway_heatmap_HC.png"),
                            width="100%", height="auto", class_name="heatmap-image"),
                        rx.heading("Pathway ErrorBar", font_size="1.2em"),
                        rx.image(
                            src=rx.asset("functional_analysis_pathway/pathway_errorbar_HC.png"),
                            width="100%", height="auto", class_name="heatmap-image"),
                        value="functional-analysis-prediction", class_name="tab-content"),
                    default_value="overview"),
                direction="column", width="100%", spacing="2", padding="1em")))
