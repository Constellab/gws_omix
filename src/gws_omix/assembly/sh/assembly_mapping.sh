# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

mapping_file=$1
chromosome_name=$2
genome_fasta=$3
filtering_threshold=$4
threads=$5

# masking uncovered region

samtools faidx -o tmp.chromosome.fasta $genome_fasta $chromosome_name ;
samtools faidx -o genome.fasta_index.fai $genome_fasta ;
samtools view -@ $threads -b -o tmp.mapping.bam -q $filtering_threshold $mapping_file ;
bedtools bamtobed -i tmp.mapping.bam | bedtools merge -i stdin | bedtools complement -L -i stdin -g <(grep "^$chromosome_name" genome.fasta_index.fai | cut -f1,2 ) | bedtools maskfasta -fi tmp.chromosome.fasta -bed stdin -fo tmp.chromosome.masked.fasta

# variant calling on the previous masked sequence (including SNPs and InDels)

samtools mpileup -uf tmp.chromosome.masked.fasta tmp.mapping.bam | bcftools call -mv -Oz -o calls.vcf.gz ;

# keep only homozygous variants

cat <(zcat calls.vcf.gz | grep "^#" )  <(zcat calls.vcf.gz | grep -P "\t1\/1:") >  calls.homozyg.vcf ; 
bcftools view calls.homozyg.vcf -Oz -o calls.homozyg.vcf.gz ;
bcftools index --csi calls.homozyg.vcf.gz ;

# file normalisation to avoid format error and strand in last step (error in ref sequence)

bcftools norm  -f tmp.chromosome.masked.fasta -m +any -Oz -o calls.norm.vcf.gz  calls.homozyg.vcf.gz ;
bcftools index --csi calls.norm.vcf.gz ;

# filtering variants according to quality

bcftools +fill-tags calls.norm.vcf.gz -Oz -o out.vcf.gz ;
bcftools filter -sLowQual -g3 -G10 -e'%QUAL<10 || (RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || %MAX(DP4)<=3 || %MAX(DP4)/%MAX(DP)<=0.3' -Oz -o test.filter.csv.gz out.vcf.gz ;
zcat test.filter.csv.gz | grep -v -P "LowQual" > exclude_low_qual.vcf ;
bcftools view exclude_low_qual.vcf -Oz -o final.vcf.gz ;
bcftools index --csi final.vcf.gz ;

# consensus sequence in fasta format

bcftools consensus -f tmp.chromosome.masked.fasta final.vcf.gz > $mapping_file"."$chromosome_name".mapping_assembly.fna" ;

mv final.vcf.gz $mapping_file"."$chromosome_name".variant_calling.vcf"

# moving output files 


#mv $mapping_file"."$chromosome_name".mapping_assembly.fna" ./assembly ;
#cp final.vcf.gz ./assembly/$mapping_file"."$chromosome_name".variant_calling.vcf" ;



