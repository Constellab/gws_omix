from gws_core import (AppConfig, AppType, ConfigParams, ConfigSpecs,
                      InputSpecs, OutputSpec, OutputSpecs, ReflexResource,
                      Task, TaskInputs, TaskOutputs, app_decorator,
                      task_decorator)
from gws_core.impl.file.file import File
from gws_core.io.io_spec import InputSpec


@app_decorator("OmixShowcaseAppConfig", app_type=AppType.REFLEX,
               human_name="Generate show case app")
class OmixShowcaseAppConfig(AppConfig):

    # retrieve the path of the app folder, relative to this file
    # the app code folder starts with a underscore to avoid being loaded when the brick is loaded
    def get_app_folder_path(self):
        return self.get_app_folder_from_relative_path(__file__, "_omix_showcase")


@task_decorator("GenerateOmixShowcase", human_name="Generate OmixShowcase app",
                style=ReflexResource.copy_style(), hide=True)
class GenerateOmixShowcase(Task):
    """
    Task that generates the OmixShowcase app.
    """

    input_specs = InputSpecs({
        'data_file': InputSpec(File, human_name="Data Excel File",
                               short_description="Excel file used for the showcase, containing data to be displayed", is_optional=False)
    })
    output_specs = OutputSpecs({
        'reflex_app': OutputSpec(ReflexResource)
    })

    config_specs = ConfigSpecs({})

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        reflex_app = ReflexResource()

        reflex_app.set_app_config(OmixShowcaseAppConfig())
        reflex_app.name = "OmixShowcase"
        reflex_app.add_resource(inputs['data_file'], create_new_resource=False)

        return {"reflex_app": reflex_app}
