#!/usr/bin/bash

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

# bash script salmon_index_cmd.sh :
#
# This script run salmon index program in salmonIndex method build_command()
# Arguments are given inside the build_command() in the cmd array
#

scriptR=$1
countTable=$2
metadataTable=$3
metadataColumn=$4
outputName=$5
threshold=$6


echo "# Performing DESEQ2 analysis : count-table: "$countTable" ; metadata-table: "$metadataTable" ; metadata-column-to-use: "$metadataColumn" ; output-file-prefix: "$outputName" ###"

Rscript --vanilla $scriptR $countTable $metadataTable $metadataColumn $outputName

for i in $outputName.*.txt ;do cat $i | awk -v var=$i -v thrld=$threshold '{ if($6<=thrld){ print $1"\t"$2"\t"var  } }' ;done | awk 'BEGIN{print "id\tfold-change\tfile_name"}{print $0}' > SummaryTable.csv

#for i in $outputName.*.txt ;do cat $i | awk -v var=$i -v thrld=$threshold '{ if($6<=thrld){ print $1"\t"$2"\t"var  } }' ;done | perl -ne 'chomp; @t=split/\t/; $cpt++; $h{$t[2]}{$t[0]}=$t[3] ; $hGenes{$t[0]}=0; $hConditions{$t[2]}=0; END{ print "id"; foreach $k (sort keys %hGenes){ print "\t"$k;} print "\n;" foreach $k (sort keys %hGenes){ foreach $l (sort keys %hConditions){ print $l; if($h{$k}{$l}){ print "\t",$h{$t[2]}{$t[0]} } else{print "\tNA"} } } }' > SummaryTable.heatmap.csv



echo "# Performing DESEQ2 analysis : count-table: "$countTable" ; metadata-table: "$metadataTable" ; metadata-column-to-use: "$metadataColumn" ; output-file-prefix: "$outputName" ; threshold-for-summary : "$threshold" ###"

echo "# DONE ###"



