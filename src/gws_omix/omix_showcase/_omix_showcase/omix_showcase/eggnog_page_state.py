import pandas as pd
import reflex as rx
from gws_core import File
from gws_reflex_main import ReflexMainState


class EggnogPageState(ReflexMainState, rx.State):

    @rx.var(cache=True)
    async def database_data(self) -> pd.DataFrame:
        resources = await self.get_resources()
        if not resources or len(resources) != 2:
            return pd.DataFrame()
        data_file: File = resources[0]
        if not data_file:
            return pd.DataFrame()
        data: pd.DataFrame = pd.read_excel(data_file.path, sheet_name=0)
        if data is None or data.empty:
            return pd.DataFrame()
        return data

    @rx.var(cache=True)
    async def constellab_data(self) -> pd.DataFrame:
        resources = await self.get_resources()
        if not resources or len(resources) != 2:
            return pd.DataFrame()
        data_file: File = resources[0]
        if not data_file:
            return pd.DataFrame()
        data: pd.DataFrame = pd.read_excel(data_file.path, sheet_name=1)
        if data is None or data.empty:
            return pd.DataFrame()
        return data

    @rx.var(cache=True)
    async def download_database_data(self) -> str:
        """Event handler to download the database data."""
        resources = await self.get_resources()
        if not resources or len(resources) != 2:
            return 'database_data.xlsx'
        data_file: File = resources[0]
        if not data_file:
            return 'database_data.xlsx'

        src_path = data_file.path
        dest_dir = rx.get_upload_dir()
        dest_path = dest_dir / 'database_data.xlsx'

        data = pd.read_excel(src_path, sheet_name=0)
        data.to_excel(dest_path, index=False)
        return 'database_data.xlsx'

    @rx.var(cache=True)
    async def download_constellab_data(self) -> str:
        """Event handler to download the Constellab data."""
        resources = await self.get_resources()
        if not resources or len(resources) != 2:
            return 'constellab_data.xlsx'
        data_file: File = resources[0]
        if not data_file:
            return 'constellab_data.xlsx'

        src_path = data_file.path
        dest_dir = rx.get_upload_dir()
        dest_path = dest_dir / 'constellab_data.xlsx'

        data = pd.read_excel(src_path, sheet_name=1)
        data.to_excel(dest_path, index=False)
        return 'constellab_data.xlsx'
