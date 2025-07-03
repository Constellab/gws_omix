#!/usr/bin/env python3
# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import subprocess
import sys


def download_and_decompress(db_name: str):
    print(f"[INFO] Downloading RefSeq database: {db_name}")
    subprocess.run(["update_blastdb.pl", db_name, "--decompress"], check=True)
    print(f"[INFO] {db_name} downloaded and decompressed.")


def move_db_files(db_name: str, output_dir: str):
    print(f"[INFO] Moving {db_name} files to {output_dir}")
    files = [f for f in os.listdir(".") if f.startswith(db_name + ".")]
    for f in files:
        os.rename(f, os.path.join(output_dir, f))
    print(f"[INFO] Files moved: {files}")


def build_refseq_db(output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    for db in ["refseq_rna", "refseq_protein"]:
        download_and_decompress(db)
        move_db_files(db, output_dir)
    print(f"[INFO] Finished. Databases are in {output_dir}")


if __name__ == "__main__":
    output_dir = sys.argv[1]
    try:
        build_refseq_db(output_dir)
    except Exception as e:
        print(f"[ERROR] {e}")
        sys.exit(1)
