import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class Translator:
    def __init__(self, input_path, fasta_output_path):
        self.input_path = input_path
        self.fasta_output_path = fasta_output_path

    def translate_sequences(self):
        try:
            # Read nucleotide sequences from input FASTA file
            sequences = SeqIO.to_dict(SeqIO.parse(self.input_path, "fasta"))

            # Translate nucleotide sequences to protein sequences
            translated_sequences = {}
            for name, seq_record in sequences.items():
                # Translate nucleotide sequence
                translated_seq = SeqRecord(seq_record.seq.translate(), id=name, description="")
                translated_sequences[name] = translated_seq

            # Write translated protein sequences to output FASTA file
            SeqIO.write(translated_sequences.values(), self.fasta_output_path, "fasta")

            print(f"Translation successful. Protein sequences saved to {self.fasta_output_path}")

        except Exception as e:
            print(f"Error: {e}")

# Example usage
#input_path = "/lab/user/data/ubiome_local_test/diamond_blast/translation_nuc_prot/P_murina.fa"
#output_fasta_path = "/lab/user/data/ubiome_local_test/diamond_blast/translation_nuc_prot/output_protein.fasta"

translator = Translator(sys.argv[1], sys.argv[2])
translator.translate_sequences()
