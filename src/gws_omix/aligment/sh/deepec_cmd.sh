#!/usr/bin/bash

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

# bash script blastp_cmd.sh :
#
# This script run DeepEC method build_command()
# Arguments are given inside the build_command() in the cmd array
#

deepec_py=$1
fasta_file=$2
working_dir=$3
output_file=$4

mkdir tmp_results ;
python $deepec_py -i $fasta_file -o ./tmp_results ;
mv ./tmp_results/log_files/DeepEC_Result_DL.txt $output_file ;

#ls ./tmp_results ;
#cat ./tmp_results/log_files/DeepEC_Result_DL.txt > $output_file ;