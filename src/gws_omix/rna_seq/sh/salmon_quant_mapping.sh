#!/usr/bin/bash

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

fastq_folder=$1
thrds=$2
metadata_file=$3
salmon_index=$4

cat $metadata_file | grep -v "^#" | sed '1d' | cut -f2,3 | awk '{print $1"#"$2;}' > tmp.metadata.tsv ;

cat tmp.metadata.tsv ;

for i in $( cat tmp.metadata.tsv )
do 
    echo $i ;
    frwd=$(echo $i | tr '#' '\t' | cut -f1);
    rvrs=$(echo $i | tr '#' '\t' | cut -f2); 
    echo $rvrs $frwd ;
    ls $fastq_folder/$rvrs $fastq_folder/$frwd ;
    salmon quant -i $salmon_index -p $thrds -l A -1 $fastq_folder/$frwd -2 $fastq_folder/$rvrs --validateMappings -o ./$frwd.quant.sf ;
    #mv quant.sf $frwd.quant.sf ;
done


salmon quantmerge --threads $thrds --quants *.quant.sf --column tpm -o salmon_quantmerge.tpm_count.txt ;
salmon quantmerge --threads $thrds --quants *.quant.sf --column numreads -o salmon_quantmerge.raw_count.txt ;
