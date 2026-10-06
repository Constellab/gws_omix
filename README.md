<p align="center">
  <img src="https://constellab.space/assets/fl-logo/constellab-logo-text-white.svg" alt="Constellab Logo" width="80%">
</p>

<br/>

# 👋 Welcome to GWS Omix

```gws_omix``` is a [Constellab](https://constellab.io) library — libraries are called bricks — developed by [Gencovery](https://gencovery.com/). GWS stands for Gencovery Web Services.

## 🚀 What is Constellab?


✨ [Gencovery](https://gencovery.com/) is a software company behind [Constellab](https://constellab.io), the leading open and secure digital infrastructure designed to consolidate data and unlock its full potential in the life sciences industry. Our mission is to give everyone access to data, so that it can improve people's health and well-being.

🌍 Constellab is free to use through our Fair Open Access offer. [Sign up here](https://constellab.space/). More information about the Open Access offer is available here (link to be defined).


## ✅ Features

`gws_omix` is dedicated to omics data analysis. It wraps widely-used bioinformatics tools as ready-to-use Constellab tasks:
- **Sequence alignment & search**: BLAST against NCBI (web) or a local RefSeq database, DIAMOND alignment with EC number mapping, multiple sequence alignment and visualisation
- **RNA-seq pipeline**: quality control (FastQC), read trimming (Trimmomatic, Fastp), genome and transcriptome indexing and mapping (STAR, HISAT2, Salmon), read counting (FeatureCounts), quality report aggregation (MultiQC) and differential expression analysis (pyDESeq2, including multi-contrast designs)
- **Functional enrichment**: over-representation analysis (ORA), gene set enrichment analysis (GSEA), GAF-to-GMT gene set conversion and gene ID conversion
- **KEGG pathway analysis**: KEGG enrichment with Pathview visualisation
- **Genome visualisation**: circular genome plots (pyCirclize) and comparative genome views (pyGenomeViz)
- **Phylogenetics**: tree construction (IQ-TREE) and tree visualisation (phyTreeViz)
- **Data acquisition**: FASTQ download from SRA and ENA
- **Showcase app**: generation of a ready-made Omix demonstration application

## 📄 Documentation

📄  The `gws_omix` brick documentation is available [here](https://constellab.community/bricks/gws_omix/latest/doc/getting-started/4ba9b3a3-aa05-4270-be76-c02f4f4d9e1e)

💫 The Constellab application documentation is available [here](https://constellab.community/bricks/gws_academy/latest/doc/getting-started/b38e4929-2e4f-469c-b47b-f9921a3d4c74)

## 🛠️ Installation

The `gws_omix` brick requires the `gws_core` brick.

### 🔥 Recommended method

The easiest way to install a brick is through the Constellab platform. Our Fair Open Access offer gives you a free cloud data lab in which bricks can be installed directly. [Sign up here](https://constellab.space/)

To learn more about the data lab, see [Overview](https://constellab.community/bricks/gws_academy/latest/doc/digital-lab/overview/294e86b4-ce9a-4c56-b34e-61c9a9a8260d) and [Data lab management](https://constellab.community/bricks/gws_academy/latest/doc/digital-lab/on-cloud-digital-lab-management/4ab03b1f-a96d-4d7a-a733-ad1edf4fb53c)

### 🔧 Manual installation

This section is for users who prefer to install the brick manually, either on their own machine or in the Constellab Codelab.

We recommend using Ubuntu 22.04 with Python 3.10.

#### Usage


▶️ Start the server:

```bash
gws server run
```

🕵️ Run a single unit test:

```bash
gws server test [TEST_FILE_NAME]
```

Replace `[TEST_FILE_NAME]` with the name of a test file (without the `.py` extension) from the `tests` folder, and run the command from the brick folder.

🕵️ Run the whole test suite:

```bash
gws server test all
```

📌 VSCode users can rely on the predefined run configurations in `.vscode/launch.json`.

## 🤗 Community

🌍 Join the Constellab community [here](https://constellab.community/) to share and explore stories, code snippets and bricks with other users.

🚩 Feel free to open an issue if you have a question or a suggestion.

☎️ You can also reach out to us through our website: [Constellab](https://constellab.io/).

## 🌎 License

```gws_omix``` is entirely free and open-source, released under the [GNU Affero General Public License v3.0](https://www.gnu.org/licenses/agpl-3.0.en.html).

<br/>


This brick is maintained with ❤️ by [Gencovery](https://gencovery.com/).

<p align="center">
  <img src="https://framerusercontent.com/images/Z4C5QHyqu5dmwnH32UEV2DoAEEo.png?scale-down-to=512" alt="Gencovery Logo"  width="30%">
</p>
