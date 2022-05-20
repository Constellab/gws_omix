# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

fasta_file=$1

bwa index $fasta_file

mv * ./bwa_index
cp $fasta_file ./bwa_index