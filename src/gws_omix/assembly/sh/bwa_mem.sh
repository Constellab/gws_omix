# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

genome_fasta_directory=$1
fastq_file_frw=$2
fastq_file_rvs=$3


bwa mem $fasta_file

mv * ./bwa_index
cp $fasta_file ./bwa_index