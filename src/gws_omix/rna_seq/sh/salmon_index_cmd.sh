#!/usr/bin/bash

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

# bash script salmon_index_cmd.sh :
#
# This script run salmon index program in salmonIndex method build_command()
# Arguments are given inside the build_command() in the cmd array
#

genomeFasta=$1
#annota=$2
threads=$2
genomeName=$3

# gffread -w transcripto.tmp.fa -g $genomeFasta $annota ;
# cat transcripto.tmp.fa  | cut -d " " -f 1 > transcripto.tmp.2.fa ;
# rm transcripto.tmp.fa ;
#grep "^>" $genomeFasta | cut -d " " -f 1 > decoys.txt ;
#sed -i.bak -e 's/>//g' decoys.txt ;
#cat transcripto.tmp.2.fa $genomeFasta > gentrome.fa.gz  ;

#mkdir salmon_index_folder

salmon index -k 31 -t $genomeFasta -p $threads -i salmon_index ; # -d decoys.txt  -t gentrome.fa.gz -

#mv *.salmon_index salmon_index_folder


