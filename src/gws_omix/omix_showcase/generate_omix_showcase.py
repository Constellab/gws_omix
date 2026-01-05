from gws_core import (
    AppConfig,
    AppType,
    ConfigParams,
    ConfigSpecs,
    InputSpecs,
    OutputSpec,
    OutputSpecs,
    ReflexResource,
    Task,
    TaskInputs,
    TaskOutputs,
    app_decorator,
    task_decorator,
)
from gws_core.config.param.param_spec import StrParam
from gws_core.impl.file.file import File
from gws_core.impl.file.folder import Folder
from gws_core.io.io_spec import InputSpec


@app_decorator("OmixShowcaseAppConfig", app_type=AppType.REFLEX,
               human_name="Generate show case app")
class OmixShowcaseAppConfig(AppConfig):

    # retrieve the path of the app folder, relative to this file
    # the app code folder starts with a underscore to avoid being loaded when the brick is loaded
    def get_app_folder_path(self):
        return self.get_app_folder_from_relative_path(__file__, "_omix_showcase")


@task_decorator("GenerateOmixShowcase", human_name="Generate OmixShowcase app",
                style=ReflexResource.copy_style())
class GenerateOmixShowcase(Task):
    """
    Task that generates the OmixShowcase app.
    """

    input_specs = InputSpecs({
        'eggnog_data_file': InputSpec(File, human_name="Eggnog gData Excel File",
                                      short_description="Excel file used for the eggnog showcase, containing data to be displayed", optional=False),
        'sixteen_s_data_file': InputSpec(Folder, human_name="16S Data Folder",
                                         short_description="Folder containing the 16S data files to be displayed", optional=False),
    })
    output_specs = OutputSpecs({
        'reflex_app': OutputSpec(ReflexResource)
    })

    config_specs = ConfigSpecs({
        'quality_check_forward_reads': StrParam(human_name='Quality Check Forward Reads URL',
                               short_description='The iframe URL for the quality check forward reads'),
        'quality_check_reverse_reads': StrParam(human_name='Quality Check Reverse Reads URL',
                               short_description='The iframe URL for the quality check reverse reads'),
        'rarefaction_analysis_observed_features': StrParam(
            human_name='Rarefaction Analysis Observed Features URL',
            short_description='The iframe URL for the rarefaction analysis observed features'),
        'rarefaction_analysis_shannon_index': StrParam(
            human_name='Rarefaction Analysis Shannon Index URL',
            short_description='The iframe URL for the rarefaction analysis shannon index'),
        'taxonomy_analysis_kingdom': StrParam(
            human_name='Taxonomy Analysis Kingdom URL',
            short_description='The iframe URL for the taxonomy analysis kingdom'),
        'taxonomy_analysis_phylum': StrParam(
            human_name='Taxonomy Analysis Phylum URL',
            short_description='The iframe URL for the taxonomy analysis phylum'),
        'taxonomy_analysis_class': StrParam(
            human_name='Taxonomy Analysis Class URL',
            short_description='The iframe URL for the taxonomy analysis class'),
        'taxonomy_analysis_order': StrParam(
            human_name='Taxonomy Analysis Order URL',
            short_description='The iframe URL for the taxonomy analysis order'),
        'taxonomy_analysis_family': StrParam(
            human_name='Taxonomy Analysis Family URL',
            short_description='The iframe URL for the taxonomy analysis family'),
        'taxonomy_analysis_genus': StrParam(
            human_name='Taxonomy Analysis Genus URL',
            short_description='The iframe URL for the taxonomy analysis genus'),
        'taxonomy_analysis_species': StrParam(
            human_name='Taxonomy Analysis Species URL',
            short_description='The iframe URL for the taxonomy analysis species'),
    })

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        reflex_app = ReflexResource()

        reflex_app.set_app_config(OmixShowcaseAppConfig())
        reflex_app.name = "OmixShowcase"
        reflex_app.add_resource(inputs['eggnog_data_file'], create_new_resource=False)
        reflex_app.add_resource(inputs['sixteen_s_data_file'], create_new_resource=False)

        reflex_app.set_params(params)

        reflex_app.set_requires_authentication(False)

        return {"reflex_app": reflex_app}
