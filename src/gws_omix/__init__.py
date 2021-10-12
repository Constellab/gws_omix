# Blast_EC
from .aligment.blast_ec import BlastEC
from .aligment.blast_get_ec_list import BlastEcGetEcList
from .file.blast_ec_file import BlastECFile



# DeepEC
from .aligment.deepec import DeepEC
from .aligment.deepec_get_ec_list import DeepEcGetEcList
from .file.deepec_file import DeepECFile
from .file.fasta_file import FastaFile



# Merge EC List
from .aligment.blast_ec_and_deepec_get_ec_list import ECListMerger
from .file.ec_list_file import ECListFile



# Gmap indexing
from .aligment.gmap_index import GmapIndex

# Gmap alignment
from .aligment.gmap_align_to_gff3 import GmapAlignGFF3
from .file.gff3_file import GFF3File
from .aligment.gmap_align_to_gtf import GmapAlignGTF
from .file.gtf_file import GTFFile
from .aligment.gmap_align_to_fasta import GmapAlignFasta
