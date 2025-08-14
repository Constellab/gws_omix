#!/usr/bin/env python3
import os
import subprocess
import sys
from contextlib import contextmanager

DBS = ["refseq_rna", "refseq_protein"]
# Any of these indicate a present/usable DB (covers protein/nucleotide + alias files)
DB_MARKERS = [".pal", ".nal", ".pin", ".nin", ".psq", ".nsq", ".pdb", ".ndb", ".phr", ".nhr"]

@contextmanager
def chdir(path):
    prev = os.getcwd()
    os.makedirs(path, exist_ok=True)
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev)

def db_exists(db_basename: str, folder: str) -> bool:
    for ext in DB_MARKERS:
        if os.path.exists(os.path.join(folder, db_basename + ext)):
            return True
    return False

def ensure_db(db_name: str, output_dir: str):
    if db_exists(db_name, output_dir):
        print(f"[INFO] {db_name} already present in {output_dir}. Skipping download.")
        return
    print(f"[INFO] Downloading {db_name} into {output_dir} ...")
    with chdir(output_dir):
        subprocess.run(["update_blastdb.pl", db_name, "--decompress"], check=True)
    if not db_exists(db_name, output_dir):
        raise RuntimeError(f"Downloaded {db_name}, but expected files not found in {output_dir}")
    print(f"[INFO] {db_name} ready in {output_dir}")

def main():
    output_dir = sys.argv[1]
    for db in DBS:
        ensure_db(db, output_dir)
    print(f"[INFO] Finished. Databases are in {output_dir}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERROR] {e}")
        sys.exit(1)
