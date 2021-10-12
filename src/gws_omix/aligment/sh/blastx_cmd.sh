#!/usr/bin/bash

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

# bash script blastp_cmd.sh :
#
# This script run blastx aligner program in blastEC method build_command()
# Arguments are given inside the build_command() in the cmd[] array
#

blastDb=$1
fastaFile=$2
eValue=$3
threadNbr=$4
coverage=$5
hitNbr=$6
outputFile=$7

blastx -db $blastDb -query $fastaFile -evalue $eValue -num_threads $threadNbr -qcov_hsp_perc $coverage -num_alignments $hitNbr -outfmt  "7 qaccver saccver pident qcovs qcovhsp length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore" -show_gis -task "blastx-fast" -out $outputFile ;
