# Dependencies

This document lists the dependencies of the `gws_omix` brick: the Constellab bricks it requires, and the external bioinformatics tools/packages used by its tasks.

## Constellab bricks

Declared in [`settings.json`](./settings.json):

| Brick | Version |
|---|---|
| `gws_core` | >= 0.22.1 |

## Python

The brick itself targets Python >= 3.8 (see [`setup.py`](./setup.py)).

Most tasks don't run in the main Python environment: each one that wraps an external bioinformatics tool creates its own **isolated environment** at runtime (via `CondaShellProxy` or `PipShellProxy`), defined by an env file (`*_env.yml` for conda, `*_env.txt` for pip) next to the task's source file. This keeps tool versions pinned and avoids conflicts between tasks (e.g. different tasks can use different Python versions). The table below lists the environments that are actually wired up to a task.

## Task environments

### Sequence alignment (`aligment/`)

| Task | Env file | Key packages |
|---|---|---|
| NCBI web BLAST (`blast_web.py`) | `base_env/blast_web_env.yml` | `biopython==1.70`, `blast==2.16.0`, `pandas` |
| RefSeqDB Local BLAST (`blast_local_task.py`) | `base_env/blast_web_env.yml` | `biopython==1.70`, `blast==2.16.0`, `pandas` |
| Build RefSeq BLAST DB (`blast_refseq_makedb.py`) | `base_env/blast_web_env.yml` | `biopython==1.70`, `blast==2.16.0`, `pandas` |
| Diamond / EC number mapping (`blast_diamond.py`, `diamond_ECnumber.py`) | none — runs in the main environment (pandas only) | — |
| Multiple Sequence Alignment (`multiple_sequence_alignment_vis.py`) | `base_env/msaviz_env.yml` | `pymsaviz==0.5.0`, `mafft==7.525`, `biopython==1.85`, `matplotlib==3.10.6` |

### RNA-seq (`rna_seq/`)

| Task | Env file | Key packages |
|---|---|---|
| FastQC (`quality_check/fastq_init.py`) | `rna_seq/quality_check/fastq_init_env.yml` | `fastqc`, `trimmomatic`, `star`, `pyyaml` |
| Trimmomatic (`trimming/trimmomatic.py`) | `rna_seq/trimming/trimmomatic_env.yml` | `trimmomatic==0.39` |
| Fastp (`trimming/fastp.py`) | `rna_seq/trimming/fastp_env.yml` | `fastp==0.24.0` |
| STAR genome index (`star_index/star_index.py`) | `rna_seq/star_index/star_index_env.yml` | `star` |
| HISAT2 index & align (`mapping_genome/hisat2_index.py`, `mapping_genome/hisat2.py`) | `rna_seq/mapping_genome/hisat2_env.yml` | `hisat2==2.2.1`, `samtools==1.21` |
| Salmon index, quant & merge matrix (`mapping_transcriptome/salmon*.py`) | `rna_seq/mapping_transcriptome/salmon_env.yml` | `salmon==1.10.3`, `pandas` |
| FeatureCounts (`featurecounts/featurecounts.py`) | `rna_seq/featurecounts/featurecounts_env.yml` | `subread==2.0.8` |
| MultiQC (`multiqc/multiqc.py`) | `rna_seq/multiqc/multiqc_env.txt` (pip) | `multiqc==1.27.1`, `numpy==2.2.4` |
| pyDESeq2 multi-contrast (`multipydesq2.py`) | `base_env/Pydesq2_env.yml` | `pydeseq2==0.5.2`, `scanpy==1.11.1`, `plotly==6.0.1` |

### Functional / GO enrichment (`go_enrichment/`)

| Task | Env file | Key packages |
|---|---|---|
| ORA analysis (`ora_analysis/ora_analysis.py`) | `go_enrichment/ora_analysis/ora_analysis_env.txt` (pip) | `gprofiler-official==1.0.0`, `pandas==2.2.2`, `matplotlib==3.9.0`, `seaborn==0.13.2` |
| GSEA (`gsea_analysis/gsea.py`) and GAF→GMT (`gsea_analysis/gaf_to_gmt.py`) | `go_enrichment/gsea_analysis/gsea_env.txt` (pip) | `gseapy==1.1.11`, `pandas==2.2.2`, `numpy==1.26.4` |
| Gene ID conversion (`file/genes_id_conversion/genes_id_conversion.py`) | `file/genes_id_conversion/genes_id_conversion_env.txt` (pip) | `gprofiler-official==1.0.0`, `pandas==2.2.2`, `requests==2.32.4` |

### KEGG (`kegg/`)

| Task | Env file | Key packages |
|---|---|---|
| KEGG enrichment + Pathview (`kegg_visualisation.py`) | `kegg/kegg_r_env.yml` | `r-base=4.4.3`, `bioconductor-pathview=1.46.0`, `pandas=3.0.0` |

### Genome visualization (`genome_visuali/`)

| Task | Env file | Key packages |
|---|---|---|
| Circular Genome Visualization (`pyCirclize.py`) | `genome_visuali/pycirclize_env.yml` | `pycirclize==1.10.1` |
| Comparative Genome Visualization (`pyGenomeViz.py`) | `genome_visuali/pygenomeviz_env.yml` | `pygenomeviz==1.6.1`, `biopython==1.80`, `mummer==3.23`, `mmseqs2==17.b804f` |

### Phylogenetics (`phylogenetic_analysis/`)

| Task | Env file | Key packages |
|---|---|---|
| Phylogenetic Tree Construction (`iqtree3/iqtree3.py`) | `phylogenetic_analysis/iqtree3/iqtree3_env.yml` | `iqtree==3.0.1` |
| Phylogenetic Tree Visualization (`phyTreeViz/phyTreeViz.py`) | `phylogenetic_analysis/phyTreeViz/phytreeviz_env.yml` | `phytreeviz==0.2.0` |

### Data acquisition (`file/`)

| Task | Env file | Key packages |
|---|---|---|
| Download FASTQ from SRA/ENA (`fastq_download_sra_ena/download_sra_ena.py`) | `file/fastq_download_sra_ena/download_sra_ena_env.yml` | `fastq-dl==4.0.1`, `pandas==2.3.3` |

## Unused legacy environments

A few environment definitions in `base_env/` are not imported by any current task (kept for reference / possible future reuse, but not part of the active dependency set): `omix_env.yml`, `Diamond_env.yml`, `subread_env.yml`, `topgo2_env.yml`, `r_env.yml`, `trimfq_env.yml`, `htseq_env_2.txt`, `pygsea_pip_env.txt`.
