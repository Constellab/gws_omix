#!/usr/bin/bash

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

# bash script gmap_gtf.sh :
#
# This script run gmap align to gtf method build_command()
# Arguments are given inside the build_command() in the cmd array
#

mapping_tool=$(samtools faidx $genomeFasta ; cat $genomeFasta".fai" | cut -f2 | awk '{res+=$0}END{if(res<4000000000){print "gmap"}if(res>=4000000000){print "gmapl"}}' )

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
mapping_tool=$(samtools faidx $genomeFasta ; cat $genomeFasta".fai" | cut -f2 | awk '{res+=$0}END{if(res<4000000000){print "gmap"}if(res>=4000000000){print "gmapl"}}' );

$mapping_tool -t $threadNbr -f 2 --npaths $hitNbr --min-identity $idt --min-trimmed-coverage $cov $crsSpecies $altStart $fullLght -D $genomeFasta".gmap_index" -d $genomeFasta".gmap_index" $fastaFile > tmp.gff3; 

gffread tmp.gff3 -T -o $outputFile ;

rm tmp.gff3;