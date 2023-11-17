---
layout: post
title: Conda Environment on Pegasus and new Acer genome
date: '2023-11-17'
categories: Coding
tags: [Coding]
---

This whole saga started when I met with Dr. Nick Kron about using perl and python scripts in Pegasus, as I am trying to adapt [Dr. Michael Studivan's bioinformatics pipeline](https://github.com/mstudiva/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt) to work on Pegasus and have been having lots of issues.

I got some great tips from Nick about how to submit jobs on Pegasus with .pl and .py scripts. There are a few important things to remember:
1. Need to change perl scripts to a different header: #! /usr/bin/env perl
2. Need to activate conda module every time you open Pegasus before submitting your job (can't activate conda during script)
3. If you download .pl scripts from online, you need to make the scripts executable (turns the text green) by running chmod +x ____.pl

I also asked Nick about which Acer genome he used for his alignment, because the Libro et al. 2013 genome Michael and I used had genes that suggest symbiont contamination within the Acer-labeled genes (i.e. photosynthetic genes). 

Nick said he used the new Vollmer lab Acer genome (Selwyn and Vollmer 2023). He also said that the main issue with aligning 3' TagSeq reads to a genome is that the UTRs are typically not present or well-annotated, so it can reduce alignment rates significantly. So, Nick recommended I use [this UTR_add_extend_GTF tool](https://github.com/danilotat/UTR_add_extend_GTF/tree/main) before running the alignment step. 

I downloaded the scripts from that GitHub link and ran them successfully on the Acer_K2_genomic.gtf file. 

Here are the genomes I downloaded and code I used:
```{bash}
## Download and format reference transcriptome

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/032/359/415/GCA_032359415.1_NEU_Acer_K2/GCA_032359415.1_NEU_Acer_K2_genomic.gtf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/032/359/415/GCA_032359415.1_NEU_Acer_K2/GCA_032359415.1_NEU_Acer_K2_genomic.fna.gz
gunzip GCA_032359415.1_NEU_Acer_K2_genomic.gtf.gz
gunzip GCA_032359415.1_NEU_Acer_K2_genomic.fna.gz

# or can run this to get everything from NCBI at once

curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_032359415.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCA_032359415.1.zip"
#then unzip it

GCA_032359415.1_NEU_Acer_K2_genomic.fna #this is the genomic FASTA file
GCA_032359415.1_NEU_Acer_K2_genomic.gtf #this is the GTF file

# I think i can directly run bowtie2 with the genomic.fna file. I just need to download the symbiont fasta files now.

# downloaded Symbiodinium, Breviolum, Cladocopium, and Durusdinium fasta files from Dr. Michael Studivan's Google Drive folder. citations are:
# Diploria host = Avila-Magana et al. 2021
# Symbiodinium = Shogushi 2018
# Breviolum = Avila-Magana et al. 2021
# Cladocopium = Camp 2022
# Durusdinium = Shogushi et al. 2021
```

The following code worked for me for the extending UTR regions:

```{bash}
## adding UTRs based on Nick Kron's recommendation for better alignment of Acer 2023 genome 

module load python/3.8.7

python3 gtf_advanced_parser.py --input ../Acer_2023/GCA_032359415.1_NEU_Acer_K2_genomic.gtf --output Acer_K2_genomic_parsed.gtf

# We finished here. Here's some info:
# Processed genes: 33794
# Processed transcripts: 28059
# Added 3'UTR onto forward strand: 13963
# Added 3'UTR onto reverse strand: 14096

# next submit the process_them.sh to bsub to extend the UTRs by a threshold of 5000 (based on what Dr. Nick Kron did)
```

But then the next step which is important is that I need to convert the updated gtf file to a fasta file, because the alignment tool Michael used (bowtie2) requires a fasta file. This is because I asked ChatGPT and it gave me this:

"Yes, you can use Bowtie2 in combination with a .gtf file to build a transcriptome index, but not directly. Bowtie2 itself is designed to build indexes from genomic DNA sequences (typically in FASTA format), and it doesn't natively support the .gtf file format, which contains gene annotation data.

However, you can use the .gtf file in conjunction with a genome sequence file to create a transcriptome index suitable for alignment with Bowtie2. The common approach is to use tools like `gffread` (part of the Cufflinks suite) or other similar utilities to first extract the transcript sequences from the genome using the .gtf file. This process generates a new FASTA file containing the transcript sequences.

Here are the general steps:

1. **Extract Transcript Sequences**: Use `gffread` or a similar tool to convert the .gtf annotations and the corresponding genome sequence (in FASTA format) into a transcriptome sequence file (also in FASTA format).

   ```bash
   gffread -w transcripts.fa -g genome.fa annotations.gtf
   ```

2. **Build Bowtie2 Index**: Once you have the transcriptome sequences in a FASTA file, you can use Bowtie2 to build an index from this file.

   ```bash
   bowtie2-build transcripts.fa transcriptome_index
   ```

This approach allows you to align RNA-Seq reads directly to the transcriptome, which can be particularly useful for tasks like expression quantification or splice site analysis."

So, I had gffread downloaded already and tried running it on the new Acer genome files, and I got this error: "no valid ID found for GFF record." That means the gtf file itself is having issues and gffread isn't recognizing the header/labels. I googled the error and found that .gtf files sometimes have 'transcript_id=""' or a blank transcript id, which is tripping up gffread because it expects something there. After googling how to fix that, I came across the tool [AGAT](https://github.com/NBISweden/AGAT/tree/master) that fixes attributes in gtf and gff files.

Now this is where I'm currently stuck: I haven't been able to run any of the AGAT scripts. They are perl scripts, but there are a lot of dependencies that I don't have installed. And that's because when I tried to install agat the "traditional" way, using bioconda, it didn't actually download any of the agat packages. I think the issue is because I was trying to use the anaconda environment/module that is on the shared program space (i.e. Pegasus creators download and maintain it). So I finally bit the bullet and downloaded anaconda3 to my nethome on Pegasus, so I could have admin privileges for updates and such. 

```{bash}
wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
bash Anaconda3-2021.05-Linux-x86_64.sh
```

Now I'm trying to install something that is supposed to make the solving dependencies process go by faster (as the "solving environment" step whenever I try to download something from Anaconda can take a looooong time).

```{bash}
conda install -n base conda-libmamba-solver
```


