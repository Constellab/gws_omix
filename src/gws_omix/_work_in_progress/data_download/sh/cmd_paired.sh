#!/usr/bin/bash

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

#Initial steps, for running qiime2 you need metadata_file and fastq_folders
## paired-end project

fastq_dir=$1
fwd=$2
rvs=$3
thrds=$4
qual=$5

max_n=$6
min_size=$7


#for i in $fastq_dir/*$fwd*.gz ;do echo $i; echo $( echo $i | sed "s/$fwd/$rvs/" ); echo -e "\n######\n" ;done


for i in $fastq_dir/*$fwd*.gz ;do echo -e "### Process" $i $( echo $i | sed "s/$fwd/$rvs/" ) "\n####\n"; rvrs=$(echo $i | sed "s/$fwd/$rvs/"); trim_galore --paired --gzip --retain_unpaired --cores $thrds --max_n $max_n --length $min_size -q $qual $i $rvrs ;done

#trim_galore --paired --gzip --retain_unpaired --cores $thrds --max_n $max_n --length $min_size -q $qual $fastq_dir/*$fwd*.gz $fastq_dir/*$rvs*.gz ;

mkdir filtered_fastq_folder
mkdir singleton_fastq_folder
#mkdir filtering_stats

mv *_val_* ./filtered_fastq_folder
mv *_unpaired_* ./singleton_fastq_folder

cat *_trimming_report* > trimming_report.txt


#trim_galore --paired --gzip --retain_unpaired --cores $thrds --max_n $max_n --length $min_size -r1 $min_size -r2 $min_size -q $qual --clip_R1 $hard_5 --clip_R2 $hard_5 --three_prime_clip_R1 $hard_3 --three_prime_clip_R2 $hard_3 $fastq_dir/*$fwd*.gz $fastq_dir/*$rvs*.gz




#for i in $fastq_dir/*$fwd*.gz ;do echo -e "trim-galore "$(basename $i )"\t"$(basename $i)"\t"$(basename $i | sed "s/$fwd/$rvs/")  ;done 
