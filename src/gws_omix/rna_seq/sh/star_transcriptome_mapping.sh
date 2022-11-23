#!/usr/bin/bash

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

fastq_folder=$1
thrds=$2
ram=$3
metadata_file=$4
genome_idx=$5
#experiment_opt=$6

mkdir star_transciptome_mapping ;

cat $metadata_file | grep -v "^#" | sed '1d' | cut -f2,3 | awk '{print $1"#"$2;}' > tmp.metadata.tsv ;

cat tmp.metadata.tsv ;

zcat $genome_idx/*.compressed.fa.gz > tmp.decompressed_genome_fasta_file.fna :
head tmp.decompressed_genome_fasta_file.fna
zcat $genome_idx/*.compressed.gtf.gz > tmp.decompressed_genome_anotation_file.gtf ;
head tmp.decompressed_genome_anotation_file.gtf

for i in $( cat tmp.metadata.tsv )
do 
    echo $i ;
    frwd=$(echo $i | tr '#' '\t' | cut -f1);
    rvrs=$(echo $i | tr '#' '\t' | cut -f2); 
    echo $rvrs $frwd ;
    ls $fastq_folder/$rvrs $fastq_folder/$frwd ;
    zcat $fastq_folder/$rvrs > tmp.rvrs.fq ;
    zcat $fastq_folder/$frwd > tmp.frwd.fq ;
    STAR --genomeDir $genome_idx --runThreadN $thrds --limitGenomeGenerateRAM $ram -outSAMtype SAM --quantMode TranscriptomeSAM --quantTranscriptomeBan Singleend --readFilesIn tmp.frwd.fq tmp.rvrs.fq --outReadsUnmapped None --twopassMode None  --outSAMstrandField intronMotif --outSAMunmapped None --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --outSAMattrRGline ID:GRPundef --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 ; # --readFilesCommand " zcat "   -outSAMtype BAM Unsorted --twopassMode Basic --outSAMunmapped Within --chimMultimapNmax 20
    samtools view -S -b Aligned.out.sam | samtools sort -@ $thrds -m 1536M -O bam -o Aligned.toTranscriptome.out.bam - ;
    rm tmp.rvrs.fq tmp.frwd.fq Aligned.out.sam ;
    mv Aligned.toTranscriptome.out.bam ./star_transciptome_mapping/$frwd".star_mapping_transcriptome.bam" ;
 
    # gffread -w transcripto.tmp.fa -g tmp.decompressed_genome_fasta_file.fna tmp.decompressed_genome_anotation_file.gtf ;
    # head transcripto.tmp.fa ;
    # cat transcripto.tmp.fa  | cut -d " " -f 1 > transcripto.tmp.2.fa ;
    # rm transcripto.tmp.fa ;

    # zcat $fastq_folder/$rvrs > tmp.rvrs.fq ;
    # zcat $fastq_folder/$frwd > tmp.frwd.fq ;

    # STAR  --runMode genomeGenerate --runThreadN  $threads --limitGenomeGenerateRAM $ram --genomeSAindexNbases 10 --genomeDir ./indexed_transcripto_folder --genomeFastaFiles tmp.decompressed_genome_fasta_file.fna --sjdbGTFfile tmp.decompressed_genome_anotation_file.gtf --sjdbOverhang $readsLength ;


    # salmon quant --threads $thrds  -t transcripto.tmp.2.fa -l A -a ./star_transciptome_mapping/$frwd".star_mapping_transcriptome.bam" -o ./star_transciptome_mapping/$frwd".star_mapping_transcriptome.count.output" ; # -l $experiment_opt

done

rm tmp.metadata.tsv ;

#salmon quantmerge --threads $thrds --quants ./star_transciptome_mapping/*.star_mapping_transcriptome.count.output --column tpm -o salmon_quantmerge.tpm_count.txt ;
#salmon quantmerge --threads $thrds --quants ./star_transciptome_mapping/*.star_mapping_transcriptome.count.output --column numreads -o salmon_quantmerge.raw_count.txt ;

echo "Done";

