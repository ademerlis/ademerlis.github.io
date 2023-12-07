---
layout: post
title: pre-filtering samples based on million reads
date: '2023-08-16'
categories: coding
tags: [coding, ch2_tempvariability]
---

I have been writing up a summary of the gene expression analysis in advance of my committee meeting, but I realized something that is inconsistent between different people's codes and I'm not sure what the baseline or "correct" way is.

I wrote up the total number of reads (and other stats) obtained for all Acer and Pcli samples:

Across both time points (day 0 and day 28 of treatment), the total number of raw reads obtained via 3’ TagSeq for A. cervicornis (N=48) and P. clivosa (N=48) were 525.6 million reads and 458.2 million reads, respectively. The average (± standard deviation) number of raw reads per sample for A. cervicornis and P. clivosa were 10.95 (± 2.67) million reads and 9.55 (± 6.88) million reads, respectively. Following trimming of low-quality base assignments, Illumina adapters, and polyA tails, reads were aligned to the A.cervicornis genome (2019 version from Baums and Kitchen) and the P. clivosa transcriptome (Avila-Magana et al. 2021). The average (± standard deviation) alignment rates for A.cervicornis and P. clivosa were 66.67% (± 6.03%) and 64.81% (± 2.90%), respectively.

So there are the raw reads, and then there is the number of aligned reads (percent_alignment x trimmed reads). 

In [Natalia's code](https://github.com/China2302/SCTLD_RRC/blob/main/05_read_counts.Rmd), she filters based on the number of reads which aligned to the genome/transcriptome, in her case a cut-off of 4 million reads. I realized in my [DESeq2 analysis of Acer](https://github.com/ademerlis/ademerlis.github.io/blob/master/_posts/2023-08-14_Ch2AcerDESeq2.md?plain=1) that I filtered out 3 samples based on this 4 million read alignment cut-off. But, in the Pcli samples, I didn't filter out any of them (I don't think), I just filtered out based on counts (see [this code](https://github.com/ademerlis/temperaturevariability2023/blob/main/gene_expression/bioinformatics/Pcli/5_Pcli_txtimport.Rmd)). But a lot more of the Pcli had low counts, including some samples that were even below 1 million counts. 

And in terms of alignment rates, can I do a Pcli cut-off at the same count level as Acer, since Acer is aligning mRNA sequences to a genome while Pcli is aligning mRNA sequences to a transcriptome (aka all the transcripts that can be made from each gene)? Does that even make sense?

Then I looked back at the SCTLD jamboree paper we did with Nikki et al ([Traylor-Knowles et al. 2021](https://www.frontiersin.org/articles/10.3389/fmars.2021.681563/full)) and from the methods it says that we filtered based on a cut-off of 2 million raw reads, not aligned reads. 

So the question is, what is that cut-off? is it based off aligned reads or raw reads?

how do we treat transcripts in this scenario (for alignment cut-offs)?
