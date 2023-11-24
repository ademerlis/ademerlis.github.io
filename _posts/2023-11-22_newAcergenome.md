---
layout: post
title: new Acer genome continued
date: '2023-11-22'
categories: Coding
tags: [Coding]
---

Previous ways I tried to download it:
```{bash}
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/032/359/415/GCA_032359415.1_NEU_Acer_K2/GCA_032359415.1_NEU_Acer_K2_genomic.gtf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/032/359/415/GCA_032359415.1_NEU_Acer_K2/GCA_032359415.1_NEU_Acer_K2_genomic.fna.gz
gunzip GCA_032359415.1_NEU_Acer_K2_genomic.gtf.gz
gunzip GCA_032359415.1_NEU_Acer_K2_genomic.fna.gz

curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_032359415.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCA_032359415.1.zip"
#then unzip it

GCA_032359415.1_NEU_Acer_K2_genomic.fna #this is the genomic FASTA file
GCA_032359415.1_NEU_Acer_K2_genomic.gtf #this is the GTF file
```

When I used curl to download the whole folder, these are the files I got:
<img width="670" alt="Screen Shot 2023-11-22 at 1 01 34 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/4d4af77d-29d1-47f7-944a-cee1348b199a">

But, in the .gtf and .gff files above, when I try to grep 'intron' or 'three_prime_UTR', I can't find any non-coding regions in these files. These are necessary for alignment. 

So, I am going to try to download them again, because when I read through the [Selwyn and Vollmer 2023 paper](https://doi.org/10.1093/g3journal/jkad232), they specifically state that the genome assembly can be accessed through NCBI BioProject: PRJNA948411 and RNA sequencing data used for annotation can be accessed from NCBI BioProject: PRJNA949884. So maybe I accidentally downloaded the wrong files and that's why they're missing introns and UTRs.

https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_032359415.1/ 

<img width="541" alt="Screen Shot 2023-11-22 at 12 59 58 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/eea43af1-0053-4bdf-b49d-9888f6901b59">

