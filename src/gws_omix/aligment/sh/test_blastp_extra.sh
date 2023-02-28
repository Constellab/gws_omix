#!/usr/bin/bash

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

# bash script blastp_cmd.sh :
#
# This script run blastp aligner program in blastEC method build_command()
# Arguments are given inside the build_command() in the cmd array
#

blastDb=$1
dirDb=$blastDb".fasta_blast_index"
fastaFile=$2
eValue=$3
threadNbr=$4
coverage=$5
hitNbr=$6
outputFile=$7

#dbPath= $blastDb"_blast_index" 
echo "####################"
echo "# inputs"
echo  "# index ID : " $blastDb".fasta" " ; # query fasta file : " $fastaFile " ; # evalue : " $eValue " ; # threads : " $threadNbr " ; # coverage : " $coverage " ; # max result hit : " $hitNbr  " ; # output file : " $outputFile " ; # index dir : " $blastDb".fasta_blast_index" " and " $(basename $blastDb)
echo "####################"
echo "####################"
echo "####################"
echo "# db files detail : "
echo "db id : "
ls  $blastDb".fasta" 
ln -s $blastDb".fasta" ./
echo "db content : "
ln -s $dirDb/* ./

echo "####################"
echo "####################"
echo "####################"
echo "Command line : " 'blastp -db $(basename $blastDb".fasta") -query $fastaFile -evalue $eValue -num_threads $threadNbr -qcov_hsp_perc $coverage -num_alignments $hitNbr -outfmt  "7 qaccver saccver pident qcovs qcovhsp length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore" -show_gis -task "blastp-fast" -out $outputFile ; '

ls $dirDb 
ls uniprot-*.fasta
#ln -s $blastDb".fasta"
#ln -s $blastDb".tab"

echo "# Running blastp to perform EC number annotation (db : " uniprot*.fasta " ) ; input fasta : " $fastaFile
blastp -db uniprot*.fasta -query $fastaFile -evalue $eValue -num_threads $threadNbr -qcov_hsp_perc $coverage -num_alignments $hitNbr -outfmt  "7 qaccver saccver pident qcovs qcovhsp length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore" -show_gis -task "blastp-fast" -out $outputFile ;
echo "# DONE : blastp "
#$threadNbr