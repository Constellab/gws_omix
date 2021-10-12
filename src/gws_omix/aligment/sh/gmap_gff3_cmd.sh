#!/usr/bin/bash

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

# bash script gmap_gff3.sh :
#
# This script run gmap align to gff3 method build_command()
# Arguments are given inside the build_command() in the cmd array
#

threadNbr=$1
hitNbr=$2
idt=$3
cov=$4
crsSpecies=$5
altStart=$6
fullLght=$7
genomeFasta=$8
fastaFile=$9
outputFile=$10
mapping_tool=$(samtools faidx $genomeFasta > tmp.txt ; cat tmp.txt | cut -f2 | awk '{res+=$0}END{if(res<4000000000){print "gmap"}if(res>=4000000000){print "gmapl"}}' ; rm tpm.txt ; );

$mapping_tool -t $threadNbr -f 2 --npaths $hitNbr --min-identity $idt --min-trimmed-coverage $cov $crsSpecies $altStart $fullLght -D $(echo $genomeFasta".gmap_index" ) -d $( echo $genomeFasta".gmap_index") $fastaFile > $outputFile ;
        