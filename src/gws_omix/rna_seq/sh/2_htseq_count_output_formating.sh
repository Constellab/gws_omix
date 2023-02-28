#!/usr/bin/bash

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

# bash script salmon_index_cmd.sh :
#
# This script run salmon index program in salmonIndex method build_command()
# Arguments are given inside the build_command() in the cmd array
#

cat $(for i in *.htseq-count.txt ;do echo $i ;done | head -1) | cut -f1 > tpm.first_column.txt

#rm $(for i in *.htseq-count.txt ;do echo $i ;done | head -1) ;

for i in *.htseq-count.txt ;do join tpm.first_column.txt $i > tmp.txt ; cat tmp.txt > tpm.first_column.txt ;done 

cat tpm.first_column.txt | grep "^__" | tr ' ' '\t' > unmapped_stats.txt
cat tpm.first_column.txt | grep -v "^__" | tr ' ' '\t'  > merged.htseq-count.txt