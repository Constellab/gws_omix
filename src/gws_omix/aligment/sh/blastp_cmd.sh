#!/usr/bin/bash

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

# bash script blastp_cmd.sh :
#
# This script run blastp aligner program in blastEC method build_command()
# Arguments are given inside the build_command() in the cmd array
#

blastDbPrefix=$1
dirDb=$blastDbPrefix".fasta_blast_index.tar.gz"
fastaFile=$2
eValue=$3
threadNbr=$4
coverage=$5
hitNbr=$6
outputFile=$7
opendataDir=$8

echo "####################"
echo "# inputs"
echo  "# index ID : " $blastDbPrefix".fasta" " ; # query fasta file : " $fastaFile " ; # evalue : " $eValue " ; # threads : " $threadNbr " ; # coverage : " $coverage " ; # max result hit : " $hitNbr  " ; # output file : " $outputFile " ; # index dir : " $blastDbPrefix".fasta_blast_index" " and " $(basename $blastDbPrefix)
echo "####################"
echo "####################"
echo "####################"
echo "# db files detail : "
echo "db id : "
ls  $blastDbPrefix".fasta.gz"  # fasta_blast_index.tar.gz tab.gz
#ln -s $blastDbPrefix".fasta.gz" ./
echo "tab file : "
ls  $blastDbPrefix".tab.gz"  # fasta_blast_index.tar.gz tab.gz
#ln -s $blastDbPrefix".tab.gz" ./
echo "db content : "
ls $dirDb 
#ln -s $dirDb/* ./

echo "####################"
echo "####################"
echo "####################"
echo "# Extract DB files : "

testFasta=$blastDbPrefix".fasta"
testTab=$blastDbPrefix".tab"
testDB=$blastDbPrefix".fasta_blast_index"

if [ -f "$testFasta" ]; then
    echo "$testFasta exists."
    ln -s $testFasta ./
else 
    echo "$testFasta does not exist."
    zcat $blastDbPrefix".fasta.gz" > $testFasta
    ln -s $testFasta ./
fi

if [ -f "$testTab" ]; then
    echo "$testTab exists."
    ln -s $testTab ./
else 
    echo "$testTab does not exist."
    zcat $blastDbPrefix".tab.gz" > $testTab
    ln -s $testTab ./
fi

if [ -d "$testDB" ]; then
    echo "$testDB exists."
    ln -s $testDB/* ./
else 
    echo "$testDB does not exist."
    tar -xzvf $dirDb --directory $opendataDir
    ln -s $testDB/* ./
fi

#ln -s $blastDbPrefix".tab"

echo "# Running blastp to perform EC number annotation (db : " uniprot*.fasta " ) ; input fasta : " $fastaFile
echo "Command line : " 'blastp -db uniprot*.fasta -query $fastaFile -evalue $eValue -num_threads $threadNbr -qcov_hsp_perc $coverage -num_alignments $hitNbr -outfmt  "7 qaccver saccver pident qcovs qcovhsp length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore" -show_gis -task "blastp-fast" -out $outputFile ; '

blastp -db uniprot*.fasta -query $fastaFile -evalue $eValue -num_threads $threadNbr -qcov_hsp_perc $coverage -num_alignments $hitNbr -outfmt  "7 qaccver saccver pident qcovs qcovhsp length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore" -show_gis -task "blastp-fast" -out $outputFile ;

echo "# DONE : blastp "

#######################################

# blastDb=$1
# fastaFile=$2
# eValue=$3
# threadNbr=$4
# coverage=$5
# hitNbr=$6
# outputFile=$7

# #dbPath= $blastDb"_blast_index" 
# for i in $blastDb"_blast_index" ;do ln -s $i ;done
# ln -s $blastDb

# blastp -db $(basename $blastDb) -query $fastaFile -evalue $eValue -num_threads 2 -qcov_hsp_perc $coverage -num_alignments $hitNbr -outfmt  "7 qaccver saccver pident qcovs qcovhsp length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore" -show_gis -task "blastp-fast" -out $outputFile ;

# #$threadNbr