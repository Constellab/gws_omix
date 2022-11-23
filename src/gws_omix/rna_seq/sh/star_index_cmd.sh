#!/usr/bin/bash

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

threads=$1
ram=$2
genomeFastaFile=$3
annotationFile=$4
readsLength=$5


mkdir indexed_genome_folder ;

zcat $genomeFastaFile > tmp.decompressed_genome_fasta_file.fna :
zcat $annotationFile | egrep '\btranscript_id\s+"[^"]+"' > tmp.decompressed_genome_anotation_file.gtf ;

STAR  --runMode genomeGenerate --runThreadN  $threads --limitGenomeGenerateRAM $ram --genomeSAindexNbases 10 --genomeDir ./indexed_genome_folder --genomeFastaFiles tmp.decompressed_genome_fasta_file.fna --sjdbGTFfile tmp.decompressed_genome_anotation_file.gtf --sjdbOverhang $readsLength ;

cp $genomeFastaFile ./indexed_genome_folder/$(basename $genomeFastaFile)".compressed.fa.gz" ;
gzip -9 -c tmp.decompressed_genome_anotation_file.gtf > ./indexed_genome_folder/$(basename $annotationFile)".compressed.gtf.gz" ;

rm tmp.decompressed_genome_fasta_file.fna tmp.decompressed_genome_anotation_file.gtf ;
