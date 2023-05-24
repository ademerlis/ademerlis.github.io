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
#~/scripts/fastqc_stresshardening2022.job
#/projects/scratch/transcriptomics/allysondemerlis/scripts/fastqc_stresshardening2022.job
#purpose: quality checking of raw RNAseq reads using FASTQC on Pegasus compute node

#BSUB -J fastqc
#BSUB -q general
#BSUB -P transcriptomics
#BSUB -o fastqc_stresshardening2022.out
#BSUB -e fastqc_stresshardening2022.err
#BSUB -n 8

and="/projects/scratch/transcriptomics/allysondemerlis" 
#this may change depending on the project

cd ${and}/projects/sequences/stresshardening2022
for SAMP in *.fastq.gz

do

module load java/1.8.0_60
${and}/programs/FastQC/fastqc \
${and}/projects/sequences/stresshardening2022/$SAMP \
--outdir ${and}/projects/stresshardening2022/fastqc_results
done


```
