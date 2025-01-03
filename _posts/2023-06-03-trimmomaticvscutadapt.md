---
layout: post
title: Trimmomatic vs. Cutadapt
date: '2023-06-03'
categories: [Bioinformatics]
tags: [Trimmomatic, cutadapt, RNA-seq]
---

Before even diving into the benefits of Trimmomatic versus Cutadapt on 3' RNA-seq data, I remembered a 2020 publication which said that adapter and low-quality base trimming was actually not even necessary before alignment, and that the "soft-clipping" removal of these short sequences during the alignment stage using Subread was more successful at not having false positives (and removing relevant data) while still removing 94% of adapter sequences. Here is the link to the paper: https://academic.oup.com/nargab/article/2/3/lqaa068/5901066

I looked up this paper on ConnectedPapers to see if any publications have come out since 2020 to either support or dispute this finding, because when I look at people's publications on RNA-seq data, it is still part of the norm to use some sort of adapter trimming software tool. 

I guess this question depends on which project I'm looking at. For the Acer CCC project, that was sequenced on NOVAseq S2 flow cell and corresponds to QuantSeq 3'mRNA-Seq cDNA libraries. For the stress-hardening project, that was done using Illumina and TagSeq. I'm not sure how different that is from the NOVAseq/QuantSeq one. 

My initial thought is to just skip the trimming step and try using Subread to align (which by the way, do I align to a transcriptome or a genome? whichever is available?). But then I could also work on trimming from the Acer_CCC samples, and compare FastQC reports for trimmed + aligned versus untrimmed + aligned. Whatever the case, right now the untrimmed samples' FastQC files are pretty shitty, so I don't even know if I can move forward with them. But I feel like I have to try, right? 

Before going further down this rabbit hole, let's circle back to the original goal, which was to determine differences between Trimmomatic and Cutadapt. From what I can tell, there aren't many differences. The only difference I can find is that Trimmomatic is newer than cutadapt (2014 vs 2011) (https://arxiv.org/pdf/2109.03625.pdf).

I guess I can run both and see which works better? 

here is the link for the available adapters on Trimmomatic: https://github.com/timflutre/trimmomatic/tree/master/adapters
