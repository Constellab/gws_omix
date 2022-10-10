# Blast_EC
from .aligment.blast_ec import BlastEC
from .aligment.blast_ec_list_extractor import BlastECListExtractor
# DeepEC
from .aligment.deepec import DeepEC
from .aligment.deepec_ec_list_extractor import DeepECListExtractor
# Merge EC List
from .aligment.ec_list_merger import ECListMerger
# Gmap alignment
#from .aligment.gmap_align_to_gff3 import GmapAlignGFF3
from .aligment.gmap_align_to_gtf import GmapAlignGTF
# Gmap indexing
from .aligment.gmap_index import GmapIndex
# data filtering
from .data_filtering.trim_galore import TrimGalore
from .file.blast_ec_file import BlastECFile
from .file.deepec_file import DeepECFile
from .file.ec_list_file import ECListFile
from .file.fasta_file import FastaFile
from .file.fastq_folder import FastqFolder
from .file.gff3_file import GFF3File
from .file.gtf_file import GTFFile
from .rna_seq.deseq2_differential_expression import DESeq2DifferentialAnalysis
from .rna_seq.salmon_index import SalmonIndex
#from .aligment.gmap_align_to_fasta import GmapAlignFasta
# Merge EC List
from .rna_seq.salmon_quant_mapping import SalmonQuantMapping
