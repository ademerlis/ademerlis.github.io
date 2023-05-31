---
layout: post
title: QC scripts for Pegasus Stress-Hardening RNA-Seq Experiment
date: '2023-05-24'
categories: coding
tags: [coding]
---

These scripts are what I want to use for QC once I get access to Pegasus and can move all my sequences onto the project space. 

```{bash}

#!/bin/bash
#~/scratch/projects/and_transcriptomics/{foldername}/scripts/fastqc/fastqc_stresshardening2022.job
#/scratch/projects/and_transcriptomics/{foldername}/scripts/fastqc/fastqc_stresshardening2022.job
#purpose: quality checking of raw RNAseq reads using FASTQC on Pegasus compute node

#BSUB -J fastqc
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o fastqc_stresshardening2022.out
#BSUB -e fastqc_stresshardening2022.err
#BSUB -n 8

and="/scratch/projects/and_transcriptomics/{foldername}/" 

cd ${and}
for SAMP in *.fastq.gz

do

module load java/1.8.0_60
module load  \
${and}/$SAMP \
--outdir ${and}/fastqc_results
done


```
