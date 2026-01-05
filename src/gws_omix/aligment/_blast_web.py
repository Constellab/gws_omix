#!/usr/bin/env python3

import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
from Bio import SeqIO

BLAST_DB = {
    "blastn": "nt",
    "blastp": "nr",
    "blastx": "nr",
    "tblastn": "nt"
}

INPUT_TYPE_REQUIRED = {
    "blastn": "nucl",
    "blastp": "prot",
    "blastx": "nucl",
    "tblastn": "prot"
}

COLUMN_NAMES_RAW = [
    "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
    "stitle", "slen", "sacc", "qstart", "qend"
]

COLUMN_NAMES_FINAL = [
    "Query ID", "Subject ID", "Accession", "Description",
    "Total Score", "Query Cover", "E value",
    "Per. Ident", "Acc. Len"
]

def sanitize_id(seq_id):
    return "".join(c if c.isalnum() or c in "_-" else "_" for c in seq_id)

def validate_params(blast_program, seq_type):
    if blast_program not in BLAST_DB:
        raise ValueError(f"Unknown BLAST program: {blast_program}")
    if seq_type != INPUT_TYPE_REQUIRED[blast_program]:
        raise ValueError(
            f"BLAST program '{blast_program}' requires input of type '{INPUT_TYPE_REQUIRED[blast_program]}', but got '{seq_type}'."
        )

def get_query_lengths(input_fasta):
    return {sanitize_id(rec.id): len(rec.seq) for rec in SeqIO.parse(input_fasta, "fasta")}

def split_fasta(input_fasta, split_dir):
    fasta_files = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_id = sanitize_id(record.id)
        out_path = os.path.join(split_dir, f"{seq_id}.fa")
        SeqIO.write(record, out_path, "fasta")
        fasta_files.append((seq_id, out_path))
    return fasta_files

def compute_query_cover(group, qlen):
    intervals = list(zip(group["qstart"].min(axis=1), group["qend"].max(axis=1)))
    intervals.sort()
    merged = []
    for start, end in intervals:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    covered = sum(end - start + 1 for start, end in merged)
    return round((covered / qlen) * 100, 2)

def run_blast(seq_tuple, output_dir, blast_program, evalue, max_target_seqs, max_hsps, pause_time, query_lengths):
    seq_id, seq_file = seq_tuple
    output_file = os.path.join(output_dir, f"{seq_id}.tsv")
    temp_output = output_file + ".tmp"

    time.sleep(pause_time)

    cmd = [
        blast_program,
        "-query", seq_file,
        "-db", BLAST_DB[blast_program],
        "-remote",
        "-evalue", str(evalue),
        "-max_target_seqs", str(max_target_seqs),
        "-max_hsps", str(max_hsps),
        "-outfmt", "6 qseqid sseqid pident length evalue bitscore stitle slen sacc qstart qend",
        "-out", temp_output
    ]

    try:
        subprocess.run(cmd, check=True)
        df = pd.read_csv(temp_output, sep="\t", names=COLUMN_NAMES_RAW)
        if df.empty:
            return

        # Total Score (comme BLAST Web)
        df["Total Score"] = df.groupby(["qseqid", "sseqid"])["bitscore"].transform("sum")

        # Query Cover (fusion HSPs comme BLAST Web)
        df["qstart"], df["qend"] = df[["qstart", "qend"]].min(axis=1), df[["qstart", "qend"]].max(axis=1)
        qlen = query_lengths.get(seq_id, 0)
        cover = 0
        if qlen > 0:
            intervals = df[["qstart", "qend"]].sort_values("qstart").values.tolist()
            merged = []
            for start, end in intervals:
                if not merged or start > merged[-1][1]:
                    merged.append([start, end])
                else:
                    merged[-1][1] = max(merged[-1][1], end)
            covered = sum(end - start + 1 for start, end in merged)
            cover = round((covered / qlen) * 100, 2)
        df["Query Cover"] = cover

        # Renommage colonnes
        df = df.rename(columns={
            "qseqid": "Query ID",
            "sseqid": "Subject ID",
            "sacc": "Accession",
            "stitle": "Description",
            "evalue": "E value",
            "pident": "Per. Ident",
            "slen": "Acc. Len"
        })

        df_final = df[[
            "Query ID", "Subject ID", "Accession", "Description",
            "Total Score", "Query Cover", "E value", "Per. Ident", "Acc. Len"
        ]]

        df_final.to_csv(output_file, sep="\t", index=False)
        os.remove(temp_output)
        print(f"{seq_id}: BLAST terminé")
    except subprocess.CalledProcessError as e:
        print(f"{seq_id}: BLAST échoué: {e}")

def main():
    input_fasta        = sys.argv[1]
    output_dir         = sys.argv[2]
    split_dir          = sys.argv[3]
    seq_type           = sys.argv[4]
    blast_program      = sys.argv[5]
    evalue             = float(sys.argv[6])
    max_target_seqs    = int(sys.argv[7])
    max_hsps           = int(sys.argv[8])
    threads            = int(sys.argv[9])
    pause_time         = float(sys.argv[10])

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(split_dir, exist_ok=True)

    validate_params(blast_program, seq_type)
    query_lengths = get_query_lengths(input_fasta)
    split_records = split_fasta(input_fasta, split_dir)

    print(f"Launching {blast_program} on {len(split_records)} sequences (remote mode)")
    with ThreadPoolExecutor(max_workers=threads) as executor:
        for seq_tuple in split_records:
            executor.submit(run_blast, seq_tuple, output_dir, blast_program, evalue,
                            max_target_seqs, max_hsps, pause_time, query_lengths)

    print(f"All BLAST results saved in: {output_dir}")

    if os.path.exists(split_dir):
        shutil.rmtree(split_dir)
        print(f"Temporary split directory deleted: {split_dir}")

if __name__ == "__main__":
    main()
