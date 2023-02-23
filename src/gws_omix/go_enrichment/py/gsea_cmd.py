import sys

# import os
import pandas as pd
from pygsea import GSEA

if len(sys.argv) > 1:
    gene_list = sys.argv[1]
    gene_universe = sys.argv[2]
    results_number = sys.argv[3]
else:
    print("No arguments were passed to the script")

# Load the gene expression data as a pandas DataFrame
g_list = pd.read_csv(gene_list, index_col=0)

# Load the gene set data as a dictionary of sets
g_universe_gmt_file = GSEA.read_gmt(gene_universe)

# Create a GSEA object and perform enrichment analysis
gsea = GSEA(g_list, g_universe_gmt_file)

with open('All_results.txt', 'w') as f:
    result = gsea.run_enrichment()
    for line in result:
        f.write(line + '\n')
    # Print the top X enriched GO terms (default, Top 10, nes: The normalized enrichment score.)
    with open('Top_list.txt', 'w') as f:
        result_sorted = sorted(result, key=lambda x: float(x.split('\t')[2]), reverse=True)
        for line in result_sorted[:10]:
            f.write(line + '\n')
