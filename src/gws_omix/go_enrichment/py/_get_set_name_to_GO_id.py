import os

from goatools import obo_parser

# Load the GO annotation database
obo_file = 'go-basic.obo'  # Replace with the path to your GO annotation file
go_db = obo_parser.GODag(obo_file)

# A list of gene set names that you want to map to GO terms
gene_set_names = ['GOBP_RIBOSOMAL_SMALL_SUBUNIT_ASSEMBLY', 'GOBP_VERY_LONG_CHAIN_FATTY_ACID_METABOLIC_PROCESS']

# Map each gene set name to its corresponding GO term
for name in gene_set_names:
    try:
        go_term = go_db.query_term(name).id
        print(f'{name}: {go_term}')
    except KeyError:
        print(f'No GO term found for {name}')

##################################
# more complete

if len(sys.argv) > 1:
    gmt_f = sys.argv[1]
    go_obo_db = sys.argv[2]
else:
    print("No arguments were passed to the script")

convert_gmt_to_go(gmt_f, go_obo_db, "db_gmt_go.txt")


def convert_gmt_to_go(gmt_file, go_obo_file, output_file):
    go_data = obo_parser.GODag(go_obo_file)
    with open(gmt_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            fields = line.strip().split('\t')
            gene_set_name = fields[0]
            gene_set_genes = fields[2:]
            go_id = None
            try:
                # Attempt to get the GO term ID for the gene set name
                go_term = go_data.query_term(gene_set_name)
                go_id = go_term.id
            except KeyError:
                # If GO term ID cannot be found, use the original gene set name
                go_id = gene_set_name
            # Write out the converted line to the output file
            f_out.write("{}\t{}\t{}\n".format(go_id, gene_set_name, '\t'.join(gene_set_genes)))
