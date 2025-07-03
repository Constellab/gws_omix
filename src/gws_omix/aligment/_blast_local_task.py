#!/usr/bin/env python3

import os
import sys
import subprocess
import pandas as pd
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor

# Mapping BLAST program to sequence type and default database basename
DB_TYPE_MAP = {
    "blastn": ("nucl", "refseq_rna"),
    "tblastn": ("prot", "refseq_rna"),
    "blastx": ("nucl", "refseq_protein"),
    "blastp": ("prot", "refseq_protein")
}

RAW_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "stitle", "sacc", "slen"
]

FINAL_COLUMNS = [
    "Query ID", "Subject ID", "Accession", "Description",
    "Total Score", "Query Cover", "E value",
    "Per. Ident", "Acc. Len"
]

def sanitize_id(seq_id):
    return "".join(c if c.isalnum() or c in "_-" else "_" for c in seq_id)

def split_fasta(input_fasta, split_dir):
    fasta_files = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_id = sanitize_id(record.id)
        out_path = os.path.join(split_dir, f"{seq_id}.fa")
        SeqIO.write(record, out_path, "fasta")
        fasta_files.append((seq_id, out_path, len(record.seq)))
    return fasta_files

def calculate_query_cover(df, qlen):
    if qlen == 0:
        return 0.0
    df["start"] = df[["qstart", "qend"]].min(axis=1)
    df["end"] = df[["qstart", "qend"]].max(axis=1)
    intervals = df[["start", "end"]].sort_values("start").values.tolist()
    merged = []
    for start, end in intervals:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    covered = sum(end - start + 1 for start, end in merged)
    return round((covered / qlen) * 100, 2)

def run_blast(seq_tuple, output_dir, blast_program, db_path, evalue, max_target_seqs, max_hsps):
    seq_id, seq_file, qlen = seq_tuple
    output_file = os.path.join(output_dir, f"{seq_id}.tsv")
    temp_output = output_file + ".tmp"

    cmd = [
        blast_program,
        "-query", seq_file,
        "-db", db_path,
        "-evalue", str(evalue),
        "-max_target_seqs", str(max_target_seqs),
        "-max_hsps", str(max_hsps),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle sacc slen",
        "-out", temp_output
    ]

    try:
        subprocess.run(cmd, check=True)
        df = pd.read_csv(temp_output, sep="\t", names=RAW_COLUMNS)
        if df.empty:
            print(f"[EMPTY] No hits for {seq_id}")
            return

        df["Total Score"] = df.groupby(["qseqid", "sseqid"])["bitscore"].transform("sum")
        df["Query Cover"] = calculate_query_cover(df, qlen)

        df = df.rename(columns={
            "qseqid": "Query ID",
            "sseqid": "Subject ID",
            "sacc": "Accession",
            "stitle": "Description",
            "evalue": "E value",
            "pident": "Per. Ident",
            "slen": "Acc. Len"
        })

        df = df[FINAL_COLUMNS]
        df.to_csv(output_file, sep="\t", index=False)
        os.remove(temp_output)
        print(f"{seq_id}: BLAST finished")
    except subprocess.CalledProcessError:
        print(f"{seq_id}: BLAST failed")

def main():
    input_fasta      = sys.argv[1]
    sequence_type    = sys.argv[2]
    blast_program    = sys.argv[3]
    db_root_dir      = sys.argv[4]
    output_dir       = sys.argv[5]
    split_dir        = sys.argv[6]
    evalue           = float(sys.argv[7])
    max_target_seqs  = int(sys.argv[8])
    max_hsps         = int(sys.argv[9])
    threads          = int(sys.argv[10])

    if blast_program not in DB_TYPE_MAP:
        raise ValueError(f"BLAST program '{blast_program}' not recognized.")

    expected_type, db_basename = DB_TYPE_MAP[blast_program]
    if sequence_type != expected_type:
        raise ValueError(f"{blast_program} expects '{expected_type}' input but got '{sequence_type}'.")

    db_path = os.path.join(os.path.abspath(db_root_dir), db_basename)
    if not any(os.path.exists(db_path + ext) for ext in [".pal", ".nal", ".pin", ".nin", ".psq", ".nsq"]):
        raise FileNotFoundError(f"Database not found at: {db_path}")

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(split_dir, exist_ok=True)

    fasta_splits = split_fasta(input_fasta, split_dir)

    with ThreadPoolExecutor(max_workers=threads) as executor:
        for seq_tuple in fasta_splits:
            executor.submit(run_blast, seq_tuple, output_dir, blast_program, db_path, evalue, max_target_seqs, max_hsps)

    print(f"All results saved in: {output_dir}")
    if os.path.exists(split_dir):
        shutil.rmtree(split_dir)

if __name__ == "__main__":
    main()
