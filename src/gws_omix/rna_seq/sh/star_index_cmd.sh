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

STAR  --runMode genomeGenerate --runThreadN  $threads --limitGenomeGenerateRAM $ram --genomeDir ./indexed_genome_folder --genomeFastaFiles  $genomeFastaFile --sjdbGTFfile $annotationFile --sjdbOverhang $readsLength

cp $genomeFastaFile ./indexed_genome_folder ;
cp $annotationFile ./indexed_genome_folder ;
