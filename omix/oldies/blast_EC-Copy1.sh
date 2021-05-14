#!/usr/bin/bash

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

blastBin=$1
blastDb=$2
fastaFile=$3
tabFile=$(echo $blastDb | sed 's/faa/tab/')
eValue=$4
threadNb=$5
coverage=$6
identity=$7
outputFile=$(echo -e $8".output")
outputFileParsed=$(echo -e $8".parsed.output")
outputFileParsedBestHit=$(echo -e $8".parsed.best_hit.output")


# bash blast_EC.sh /Users/deoliveirar/anaconda3/bin/blastp /Users/deoliveirar/data/uniprot-taxonomy_4751.fungi.fasta

$blastBin -db $blastDb -query $fastaFile -evalue $eValue -num_threads $threadNb -qcov_hsp_perc $coverage -num_alignments 10 -outfmt  "7 qaccver saccver pident qcovs qcovhsp length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore" -show_gis -task "blastp-fast" -out tmp.outputfile.csv 

cat tmp.outputfile.csv | tr '\|' '#' > tmp.outputfile.2.csv

rm tmp.outputfile.csv

# Adding EC number(s) information

perl getEcNbForBlast.pl $tabFile tmp.outputfile.2.csv > $outputFile 

rm tmp.outputfile.2.csv

# Producing parsed files: (1) filtered by coverage and identity/similarity, (2) keeping only the best hit

perl blastOutputParsing.pl $identity $outputFile  > $outputFileParsed

cat $outputFileParsed | perl blastOutputParsingBestHit.pl -  | egrep "^1\t" | cut -f2- > $outputFileParsedBestHit

