#!/usr/bin/bash

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

# bash script blastp_cmd.sh :
#
# This script run DeepEC method build_command()
# Arguments are given inside the build_command() in the cmd array
#

fastaFile=$1
workingDir=$2
outputFile=$3

mkdir tmp_results ;

python deepec.py -i $fastaFile -o ./tmp_results ;

ls ./tmp_results ;

cat ./tmp_results/log_files/DeepEC_Result_DL.txt > $outputFile ;