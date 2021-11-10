
import os

from gws_core import Requests, Settings, SysProc, Zip


class UniprotDBHelper:

    #BASE_SOURCE_URL = "https://storage.gra.cloud.ovh.net/v1/AUTH_a0286631d7b24afba3f3cdebed2992aa/opendata"
    TAXONOMY = {
        "viruses": 10239,
        "bacteria": 2,
        "archaea": 2157,
        "fungi": 4751,
        "viridiplantae": 33090,
        "metazoa": 33208,
        "chordata": 7711,
        "mammalia": 40674,
        "eukaryota": 2759,
    }
    DESTINATION_DIR = "/data/gws_omix/uniprot/uniprot_kb/"

    @classmethod
    def get_taxonomy_file_path(taxonomy: str):
        return os.path.join(
            cls.DESTINATION_DIR,
            taxonomy,
            cls.get_taxonomy_file_name(taxonomy)
        )

    @classmethod
    def get_taxonomy_file_name(taxonomy: str):
        return taxonomy+".uniprotKB.faa"

    @classmethod
    def upload_db(cls, taxonomy):
        # ToDo

        settings = Settings.retrieve()
        BASE_SOURCE_URL = settings.get_variables("gws_omix:uniprot_db_url")
        source_url = BASE_SOURCE_URL + f"{taxonomy}.zip"

        cmd[
            "curl", source_url, " | ", "unzip "
        ]

        proc = SysProc.popen(
            cmd,
            cwd=cls.DESTINATION_DIR,
            stdout=subprocess.PIPE
        )

        if not os.path.exists(cls.DESTINATION_DIR):
            os.makedirs(cls.DESTINATION_DIR)

        dest_filename = ""
        Requests.download(
            url=source_url,
            dest_dir=cls.DESTINATION_DIR,
            dest_filename=taxonomy
        )

        Zip.unzip(zipfile_path)

        pass
