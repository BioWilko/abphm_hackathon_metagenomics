# Simulated Metagenome Workflow

This repository contains a Nextflow workflow for simulating metagenomic data based on reference genomes. The workflow consists of two processes, `fetch_refs` and `gen_reads`, which fetch reference genomes using their accession numbers and generate simulated reads from the reference genomes, respectively.

This work has been developed by a group of participants at the 9th Microbes and Food Safety Bioinformatics Hackathon held in Cambridge, UK in May 2023.

## Requirements

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

## Usage

1. Clone this repository:

```bash
git clone https://github.com/BioWilko/simulated-metagenome-workflow.git
cd simulated_metagenome_pipe
```

## Code Explanation

**fetch_refs process**: This process fetches reference genomes using their accession numbers. The input is a row from the metadata manifest file, and the output is a tuple containing the row, the reference genome FASTA file, and stdout. The `ref_fetch.py` script is called with the accession number, email, total reads, and proportion as arguments.

**gen_reads process**: This process generates simulated reads from the reference genomes using the `badread` package. The input is a tuple containing the row from the metadata manifest file, the reference genome FASTA file, and the number of reads to generate. The output includes a tuple containing the row, the reference genome FASTA file, and taxonomic metadata, as well as a FASTQ file with the simulated reads for each taxon.

**workflow**: This part of the code connects the `fetch_refs` and `gen_reads` processes using the `manifest_ch` channel. The output FASTQ files from the `gen_reads` process are collected into a single file called "simulated_metagenome.fastq" and stored in the specified output directory.

2. Run the Nextflow Pipeline. The following is an example
```bash
nextflow run /file/path/to/simulated-metagenome-workflow/main.nf --meta_manifest /file/path/to/manifest/file --out_dir /file/path/to/out_dir 
```
Optional arguments include `--email`, `--error_model`, `--read_len`, and `total_reads`
