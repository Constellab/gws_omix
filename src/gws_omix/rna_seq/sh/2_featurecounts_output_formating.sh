# bash script salmon_index_cmd.sh :
#
# This script run salmon index program in salmonIndex method build_command()
# Arguments are given inside the build_command() in the cmd array
#

#cut -f1,7,8,9,10,11,12 featurecounts-output.txt > featurecounts.gene_expression_matrix.txt

cut -f1,7- featurecounts-output.txt | grep -v "^#" > featurecounts.gene_expression_matrix.txt

#cat featurecounts-output.txt.summary > featurecounts.summary.txt

