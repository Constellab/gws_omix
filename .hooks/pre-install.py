

# import os

# from utils._requests import Requests
# from utils._settings import Settings
# from utils._zip import Zip

# # **********************************************************
# #
# # /!\ DO NOT ALTER THIS SECTION
# #
# # **********************************************************

# __cdir__ = os.path.dirname(os.path.realpath(__file__))
# DATA_PATH = "/data"


# def current_brick_path() -> str:
#     """ Get the current brick path """
#     return os.path.abspath(os.path.join(__cdir__, '../../'))


# def current_brick_name() -> str:
#     """ Get the current brick name """
#     return read_settings()["name"]


# def read_settings() -> dict:
#     """ Read the settings file of the current brick """
#     return Settings.read()


# def unzip(zipfile_path, output_path: str = None) -> str:
#     """ Unzip a file and return the destination path """
#     return Zip.unzip(zipfile_path, output_path)


# def download(url: str, dest_path: str) -> str:
#     """ Download a file from a remote url """
#     return Requests.download(url, dest_path)


# # **********************************************************
# #
# # UPDATE FUCNTION call_hook() TO ADD YOUR HOOK
# #
# # **********************************************************
# def call_hook():
#     """ Call hook """

#     settings = read_settings()

#     # Pull BLASTdb
#     # UniProtKB Fungi

#     url = settings["variables"]["gws_omix:uniprot_fungi_fasta_url"]
#     dest_path = settings["variables"]["gws_omix:uniprot_fungi_fasta_file"]
#     download(url, dest_path)

#     url = settings["variables"]["gws_omix:uniprot_fungi_blastdb_url"]
#     dest_path = settings["variables"]["gws_omix:uniprot_fungi_blastdb_file"]
#     download(url, dest_path)

#     url = settings["variables"]["gws_omix:uniprot_fungi_tab_url"]
#     dest_path = settings["variables"]["gws_omix:uniprot_fungi_tab_file"]
#     download(url, dest_path)


# if __name__ == "__main__":
#     call_hook()
