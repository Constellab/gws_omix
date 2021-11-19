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
idt=$2
cov=$3
hitNbr=$4
crsSpecies=$( echo $5 | awk '{if($0 == "Y"){ print "--cross-species" }else{ print "  " }}')
altStart=$( echo $6 | awk '{if($0 == "Y"){ print "--alt-start-codons" }else{ print "  " }}')
fullLght=$( echo $7 | awk '{if($0 == "Y"){ print "--fulllength" }else{ print "  " }}')
genomeFasta=$8
genome_dir=$(echo -e $genomeFasta".gmap_index")
genome_index_name=$(echo -e $( basename $genomeFasta)".gmap_index")
fastaFile=$9
outputFile=${10}
#genome_dir=${11}
#genome_index_name=$(echo -e $( basename ${11})".gmap_index")
#genome_index_name=${12}

echo $outputFile $fastaFile


mapping_tool=$(samtools faidx $genomeFasta > tmp.txt ; cat tmp.txt | cut -f2 | awk '{res+=$0}END{if(res<4000000000){print " gmap "}if(res>=4000000000){print " gmapl "}}' ; rm tmp.txt ; );


$mapping_tool $crsSpecies $altStart $fullLght -t $threadNbr -f 2 --npaths $hitNbr --min-identity $idt --min-trimmed-coverage $cov -D $genome_dir -d $genome_index_name $fastaFile | grep -v "^# Generated by GMAP" > $outputFile ;
        