#!/usr/bin/bash

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

# bash script salmon_index_cmd.sh :
#
# This script run salmon index program in salmonIndex method build_command()
# Arguments are given inside the build_command() in the cmd array
#

bam_folder=$1
annotation=$2
thrds=$3

htseq-count --help

for i in $bam_folder/*.bam
do
    htseq-count --stranded=no -m intersection-nonempty $i $annotation | awk -v sample=$(basename $i) 'BEGIN{print "Name\t"sample}{print $0;}'  > $(basename $i).htseq-count.txt ; 
    
done
