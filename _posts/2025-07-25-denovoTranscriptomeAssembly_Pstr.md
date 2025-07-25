---
layout: post
title: denovo Transcriptome Assembly for Pstr
date: '2025-07-25'
categories: Analysis, Processing
tags: [de novo transcriptome, Trinity, Pseudodiploria strigosa, Chapter 3 Reciprocal Transplant]
---

I don't have a reference genome or transcriptome to align the Chapter 3 reads of *P.strigosa* to for the reciprocal transplant study, so I need to generate a *de novo* transcriptome using the reads from the study.

However, I can't find a consensus on how many samples from the study to use. I have 192 sequence files - a forward and reverse read for each sample, or which there is 96 of them. 

Should I be using all of them? Won't that be time/resource intensive? Is that even necessary?

## GitHub Workflows
There are a few GitHub workflows I have been looking at for using Trinity:
1) [Dr. Jill Ashey](https://github.com/hputnam/Apulchra_genome/blob/main/RNA_Seq_Info/2023-08-31-Acropora-pulchra-denovo-transcriptome.md)
2) [Roberts lab](https://robertslab.github.io/resources/bio_Transcriptome-assembly/)
3) [Trinity GitHub Repository](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Assembly-Quality-Assessment)
4) [Trinity de novo Transcriptome assembly workshop](https://github.com/trinityrnaseq/RNASeq_Trinity_Tuxedo_Workshop/wiki/Trinity-De-novo-Transcriptome-Assembly-Workshop)
5) [Shrimp project](https://github.com/MaineINBRE/Trinity2.8.4Marconi/blob/master/assembly.md)

None of the above workflows specify if there is a min/max amount of sequence files to use. Dr. Jill Ashey used 12 samples (6 forward, 6 reverse). 

From [Grabherr et al. 2011 (the first Trinity publication)](https://www.nature.com/articles/nbt.1883#MOESM13), they state that 50 M pairs of reads is enough for Trinity to fully reconstruct 86% of annotated transcripts. The authors say that after 50 M paired reads, it is "saturated" or doesn't require any more data. 

[Haas et al. 2014 (the newer Trinity publication)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3875132/) doesn't talk about this. 

## Coral papers that construct *de novo* transcriptomes




