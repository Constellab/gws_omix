#!/usr/bin/bash

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

fastq_folder=$1
thrds=$2
ram=$3
metadata_file=$4
genome_idx=$5

mkdir star_transciptome_mapping ;

for i in $( cat $metadata_file | grep -v "^#" | grep -v "^sample-id\t" | cut -f2,3 | awk '{print $1"#"$2;}') ;do echo $i; frwd=$(echo $i | tr '#' '\t' | cut -f1); rvrs=$(echo $i | tr '#' '\t' | cut -f2) ; STAR --genomeDir $genome_idx --runThreadN $thrds --limitGenomeGenerateRAM $ram  -outSAMtype BAM Unsorted --quantMode TranscriptomeSAM --quantTranscriptomeBan Singleend --readFilesIn $fastq_folder/$frwd $fastq_folder/$rvrs --outReadsUnmapped None --twopassMode Basic --readFilesCommand "gunzip -c" --outSAMstrandField intronMotif --outSAMunmapped Within --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --outSAMattrRGline ID:GRPundef --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 ; rm Aligned.out.sam ; mv Aligned.toTranscriptome.out.bam ./star_transciptome_mapping/$frwd".star_mapping_transcriptome.bam" ;done

echo "Done";