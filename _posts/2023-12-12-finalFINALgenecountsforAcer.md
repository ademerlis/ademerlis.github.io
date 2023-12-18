---
layout: post
title: Acer gene counts using samtools
date: '2023-12-07'
categories: Coding
tags:
  - Coding
published: true
---

## Finally getting the gene counts for Acer samples from Chapter 2

Following the last section of [Michael Studivan's code](https://github.com/mstudiva/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt). 

```{bash}
# NOTE: Must have a tab-delimited file giving correspondence between contigs in the transcriptome fasta file and genes

cat Acer/Locatelli_2023/Acervicornis_seq2iso.tab Symbiodinium/Mstudiva_Symbiodinium_annotated_transcriptome_code/Symbiodinium_seq2iso.tab > Host_concat_seq2iso.tab
```

Job: use sam tools and [samcount.pl](https://github.com/mstudiva/tag-based_RNAseq/blob/master/samcount.pl).

arg1: SAM file (by cluster, contig, or isotig)
arg2: a table in the form 'reference_seq<tab>gene_ID', giving the correspondence of 
reference sequences to genes. With 454-deived transcriptome, the gene_ID would be isogroup; 
with Trinity-derived transcriptiome,it would be component.
  
paths:
  1. /scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer/bowtie2align_LocatelliShoguchi/sam_files/*.sam
  2. /scratch/projects/and_transcriptomics/genomes/Host_concat_seq2iso.tab
  
I need to write a for-loop to each .sam file to do this: print "samcount.pl $file $tab aligner=bowtie2 >$file.counts" but written in bash not perl.
  
```{bash}
#! /usr/bin/env bash

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir=

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer/bowtie2align_LocatelliShoguchi/sam_files"

data=($(ls *.sam))

for samp in "${data[@]}" ; do \

#build script
echo "making sam_counts script for ${samp}..."
echo "
#! /usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${samp}_samcounts
#BSUB -e ${and}/Ch2_temperaturevariability2023/3_genecounts/logs/${samp}_samcounts.err
#BSUB -o ${and}/Ch2_temperaturevariability2023/3_genecounts/logs/${samp}_samcounts.out
#BSUB -W 12:00
#BSUB -n 8
#BSUB -q general

cd \"/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer/bowtie2align_LocatelliShoguchi/sam_files\"

module load samtools/1.3

perl samcount.pl ${samp} /scratch/projects/and_transcriptomics/genomes/Host_concat_seq2iso.tab aligner=bowtie2 >${samp}.counts

" > ${and}/Ch2_temperaturevariability2023/3_genecounts/${samp}_samcounts.job

bsub < ${and}/Ch2_temperaturevariability2023/3_genecounts/${samp}_samcounts.job

done
```

When I ran this though on Pegasus, the counts files came out with nothing in them. Moreover, all the .err files have errors like:

<img width="959" alt="Screen Shot 2023-12-15 at 11 33 09 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/af05892b-22cd-4265-8d3f-d37aa7cd5e33">

Why didn't any of this work? It looks like it maybe has to do with the isogroup assignment? 

This is what Michael ran on his samples before making the seq2iso.tab file: ([link](https://github.com/mstudiva/Acropora-cervicornis-annotated-transcriptome/blob/main/tagSeq_TranscriptomeAnnotation_README.txt))

```{bash}
# Renaming gene identifiers for ease
sed -i 's/comp/Acropora/g' Acervicornis.fasta
sed -i 's/isogroup/Acropora/g' Acervicornis.fasta
```

I didn't initially run this on my fasta file but I'll try it now and see if it helps.

```{bash}
# Renaming gene identifiers for ease
sed -i 's/comp/Acropora/g' Acer2023.fasta > Acer2023_edited.fasta
sed -i 's/isogroup/Acropora/g' Acervicornis.fasta
```








