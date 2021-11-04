import sys
import json
import os
import re
import csv


def create_filtered_output_file(blast_output_file, tabular_file, id):
    gene_ec={}
    hit_parsed={}
    filtered_fp={}

    with open(tabular_file, 'r') as lines: # Create dict. containing genes with their corresponding EC number(s)
        li=lines.readlines()
        for line in li:
            if re.match("^#",line):
                pass
            else:
                li_split=line.split("\t")
                if re.match("^$",str(li_split[7])):
                    gene_ec[li_split[0]]="NA\n"                   
                else:
                    gene_ec[li_split[0]]=li_split[7]                    

    

    with open(blast_output_file, 'r') as raw_fp: 
        # Create dict. containing for each lines of the blast output 
        # (which are over the identity threshold): Hit gene's EC numbers and Best hit information
        li=raw_fp.readlines()
        best_hit_lines={}
        cpt=0

        for line in li:
            if re.match("^#",line):
                pass
            else:
                li_split=line.split("\t")

                if float(li_split[2]) >=   float(id) : # Parsing blast hit according to the identity threshold
                    cpt+=1
                    hit_gene_ids = li_split[1]
                    gene_uniprotKB_ID = hit_gene_ids.split('|')
                    gene_name = str(gene_uniprotKB_ID[1])
                    key = str(li_split[0])

                    if  key in best_hit_lines : # Give information about the best hit for each assessed gene -> Output dict.
                        hit_parsed = "{}\t{}\t{}".format(line.rstrip(),"SECONDARY_HITS",gene_ec[gene_name])
#                        print(hit_parsed)
                    else:
                        best_hit_lines[key] = 1
                        hit_parsed = "{}\t{}\t{}".format(line.rstrip(),"BEST_HIT",gene_ec[gene_name])
#                        print(hit_parsed)                                           
                    filtered_fp[cpt]=hit_parsed
                else:
                    pass

    return filtered_fp


res=create_filtered_output_file(sys.argv[1], sys.argv[2], sys.argv[3])
for rr in res:
    print(res[rr].strip())