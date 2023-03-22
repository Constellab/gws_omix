import sys

import gseapy as gp
import pandas as pd

if len(sys.argv) > 1:
    gene_list_file = sys.argv[1]
    gene_universe_file = sys.argv[2]
    results_number = sys.argv[3]
    thrds = int(sys.argv[4])
else:
    print("No arguments were passed to the script")

# Load the gene list from the gene list file
gene_list = pd.read_csv(gene_list_file, header=None)
gene_list = gene_list[0].tolist()

# Load the gene sets from the GMT file
gene_sets = gp.gsea.read_gmt(gene_universe_file)

# gene_sets = pd.read_csv(gene_universe_file, header=None)
# gene_sets = gene_universe_file

# Perform the enrichment analysis
results = gp.enrichr(gene_list=gene_list,
                     gene_sets=gene_sets,
                     organism='custom',
                     cutoff=0.5,
                     processes=thrds)

# Write the full results to a file
results.res2d.to_csv('All_results.txt', sep='\t')

# Write the top results to a separate file
top_results = results.res2d.head(int(results_number))
top_results.to_csv('Top_list.txt', sep='\t')
