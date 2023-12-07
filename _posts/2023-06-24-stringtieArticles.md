---
layout: post
title: stringtie articles
date: '2023-06-24'
categories: coding
tags: [coding]
---

I just wanted to put all these somewhere for when I need them. Stringtie is what Jill Ashey used following STAR alignment. In the past I used Subread (I think?) for Pdam. 

I think Natalia just downloaded the STAR results and used that: https://github.com/China2302/SCTLD_RRC/blob/main/06a_gene_expression_analysis.Rmd

Anyways here are the stringtie articles I found to look at later:
1. https://ccb.jhu.edu/software/stringtie/
2. https://www.nature.com/articles/nbt.3122
3. https://github.com/gpertea/stringtie
4. https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1910-1
5. https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08465-0 (review of transcriptome analysis methods with reference genome)

From Liu et al. 2022 (review of transcriptome analysis methods), they said that tools like Stringtie and Cufflinks are good for improving mapping rates because: "In some cases, the original reads may be spliced and associated with software-constructed transcriptomes to improve the alignment. The tools used for these procedures, including StringTie [8] and Cufflinks [9], can detect de novo transcripts. Moreover, when the annotation for the reference genome is incomplete, these tools can effectively fill the gap for the missing annotation information."

Interestingly, tools like Kallisto and Salmon do "pseudo-alignment":  "Several research teams recently have introduced pseudo-alignment or “alignment-free” tools. These tools, including Kallisto [10] (Fig. 1a) and Salmon [11], can directly associate the raw sequencing reads with the transcript and evaluate the gene or transcript expression levels."

But, Young et al. 2020 did STAR alignment and then Salmon for quantification: https://github.com/benyoung93/innate_immune_response_acropora_palmata/blob/master/SALMON_quant_apal.sh

So I guess you can do both? 

Here is a good summary of quantification tools from Liu et al 2022:

"Previous studies have shown that quantification tools have a greater impact on the final DE results than alignment tools [12, 13]. Commonly used quantification tools include Rcount [14], HTseq [15], StringTie [8], and Cufflinks [9] (Fig. 1a). These tools can be divided into two groups according to the evaluation standards for gene expression which can be based on counts or fragments per kilobase of transcript per million mapped reads (FPKM) values. Rcount, HTseq, and Kallisto are based on counts, while StringTie and Cufflinks are based on FPKM values. Both HTseq and Rcount count the reads mapped unambiguously to a single gene. HTseq discards the reads aligned to multiple positions and those that overlap with more than one gene [15], while Rcount assigns weights to each alignment of a multiread [14]. Therefore, Rcount is better at counting multireads and gene overlapping regions. Generally, when the reads in a dataset have good quality and length, unambiguous reads account for the majority of the transcriptome. StringTie and Cufflinks were developed by the same research team [8, 9]. Both quantify the gene expression levels based on FPKM values. The expression values for different transcripts can be determined from the results of these two programs. The resulting values for different transcripts of the same gene can be combined to obtain the gene expression values for DE analysis."

I think I'll just try using Stringtie because that's what Jill used with her Acer samples. 

