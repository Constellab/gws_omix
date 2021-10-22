#!/usr/bin/bash

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

# bash script blastp_cmd.sh :
#
# This script run blastp aligner program in blastEC method build_command()
# Arguments are given inside the build_command() in the cmd array
#

fastqFile1=$1
fastqFile2=$2
threadNbr=$3
phredThreshold=$4
unknownNucl=$5
pairedEnd=$( echo $6 | awk '{if($0 == "YES"){ print " --paired " }else{ print "  " }}')
outputDirectory=$7

# --gzip --paired --quality --output_dir

trim_galore --cores $threadNbr --quality $phredThreshold --output_dir $outputDirectory --max_n $unknownNucl $pairedEnd $fastqFile1 $fastqFile2
