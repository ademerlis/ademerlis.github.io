---
layout: post
title: PCR duplicates
date: '2023-08-21'
categories: coding
tags: [coding]
---

I met with Dr. Michael Studivan last week to discuss plans for DNA and RNA extractions of the Ch 3 Pstr samples. During this meeting, we talked about the results from the Ch 2 RNA-seq data, and he was concerned with how high the number of passed reads there were post-filtering step (see [the multiqc report](https://github.com/ademerlis/temperaturevariability2023/tree/main/gene_expression/bioinformatics/QC#2-multiqc-reports-of-trimmed-reads)). He said that when he has run TagSeq analysis, he usually get around 60% loss of reads when going from raw reads to trimmed reads. Whereas for mine, I got >95% retention of reads which passed the fastp filter. So what's happening here? Did fastp not work? 

Michael said that TagSeq is know to have a lot of PCR duplicates, and that this is important to filter out so there is only one copy per transcript.

Maybe that is what is happening here in the "sequence duplication levels" figure from multiqc? 

This first graph is from the raw reads:

<img width="816" alt="Screen Shot 2023-08-21 at 10 44 36 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/ec1a0d48-633e-4e48-aa90-867d70e441ba">

There is no graph like that for the fastp-trimmed reads. But, in the first table of general statistics, there is a column for percent duplication. I sorted it by highest to lowest, and it looks like all the Pcli samples have > 60% duplication. But, that doesn't seem to have been trimmed, because the %PF (passed filter) is ~97%. 

<img width="826" alt="Screen Shot 2023-08-21 at 10 47 14 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/3909acc9-bbde-47cf-ac5f-814441476857">

For Acer, the duplication rates are "lower" (30-50%), but the %PF is still really high.

<img width="810" alt="Screen Shot 2023-08-21 at 10 48 24 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/1bd918db-093c-4fd9-b035-7a9fa53c9c9e">

Going back through my [Ch4_AcerCCC multiqc reports](https://github.com/ademerlis/AcerCCC#2-trimmed-reads-adapters-script), it looks like TrimGalore/cutadapt actually trimmed a good chunk of base pairs (28-80%). 

![image](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/edb12491-88a1-42d0-83d8-3936bb3a8158)

But, there are still high duplication levels I think post-trimming.

![image](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/9d1f0aad-dccb-4275-ae5c-251f6cf6f93f)


**12pm meeting with Michael**:

- found several different genomes/transcriptomes for Acer, Pcli, symbiodinium, breviolum, then ran alignment tests to all the ones out there
- alignment tests:
   - first row is overall alignment, second row is percentage of reads that aligned multiple times
   - he also used the mRNA fasta files from the Baums Acer and Apal, not the original fasta file that I used.
   - so low alignment rates from Baums could be an artefact of the fact that our samples have everything (more than just mRNA), but you're trying align to a file of just mRNA.
- majority of alignments for Acer were for Symbiodinium, and majority for Pcli were Breviolum. So, that is what he used for downstream gene expression analyses. But, you wouldn't use those percent alignment rates as reporting quantitative data for symbiont typing (still need to do qPCR).
- so the Baums only had 700 genes that aligned and had counts for Acer. That's really low and suspect, so them he tried to align Polato genome. But only 159 genes that mapped. Libro alignment had ~50,000 genes. So that's good.
- for symbiont alignments, they are all going to be low percentages, and that's normal. The best one was the Shogushi genome.
- Bowtie2 -- can feed 2 different transcriptome assemblies, and it will pick the best alignment.
- for symbiodinium, got 41,000 genes.
- overall feel most confident in the Libro for Acer.
- For Pcli it's a different story:
    - aligned to Avila-Magna 2021 very high alignment rates, but lost a lot initially in trimming.
    - 50,000 genes for host, and 28,000 genes for symbiont.
- for Ch 3 -- read Kenkel et al. 2017 Past RT exp paper. Use DAPC analysis, perfect for RT. it is a way to statistically test whether transplantation resulted them being more similar to the opposite group (home:home vs. home:away).


**August 23 updates**:

I don't think the duplication plots from FastQC that I attached in the above sections are very meaningful. I came across some good explanations on the UC Davis Genome Center website that have been helping me better understand everything. Here is [their explanation on why FastQC duplication isn't important to pay attention to](https://dnatech.genomecenter.ucdavis.edu/faqs/why-does-fastqc-show-unexpectedly-high-sequence-duplication-levels-pcr-duplicates/). They specifically say: "FASTQC only analyses 50 nt of the first 100,000 reads for each file for the duplication analysis and extrapolates the dedication rates from this limited number of reads." 

And I believe that is different than the issue of PCR duplication that we are concerned with here. 

But, what is concerning and difficult to understand is the concept of PCR duplication in RNA-seq data. While Michael said it is important to remove PCR duplicates from RNA-seq data, because it is artificially enhancing the amount of counts for random reads, the UC Davis Genome Center [article says no you shouldn't remove PCR duplicates](https://dnatech.genomecenter.ucdavis.edu/faqs/should-i-remove-pcr-duplicates-from-my-rna-seq-data/) unless you do RAD-seq. But, they say that if there are unique molecular identifiers (UMIs), it is possible to remove PCR duplicates and may be helpful in some scenarios. 

I am not sure if I have UMIs in my data, but Michael refers to the Illumina adapters and headers a lot. Those could be similar concepts (things that are added during cDNA library prep to help distinguish between different samples). 

[This article from molecularecologist.com](https://www.molecularecologist.com/2016/08/25/the-trouble-with-pcr-duplicates/) linked by UC Davis Genome Center does a good job of explaining how PCR duplicates occur (but in the context of RAD data, although I think the process is the same for RNA-seq):

"There are two main steps in RAD prep that are important in generating PCR duplicates: fragmentation and PCR amplification. During fragmentation DNA is sheared randomly. PCR is then performed on these fragments to incorporate the correct adaptors for sequencing and to generate sufficient library quantity. Typical library prep uses around 10 PCR cycles. This means that the majority of your library consists of reads generated via PCR. That is, basically all of the fragments are identical to at least one other fragment in the library. This means that we will almost certainly sequence some fragments that came from the same starting DNA. Moreover, the more cycles of PCR, the higher percentage of PCR duplicates you’ll likely see."

The conclusion from the UC Davis Genome Center is that "Several studies (among them Parekh et al.  2016; below) have shown that retaining PCR- and Illumina clustering duplicates does not cause significant artifacts as long as the library complexity is sufficient.  The library complexity is in most cases directly related to the amount of starting material available for the library preparation. Chemical inhibitors present in the sample could also cause low conversion efficiency and thus reduced library complexities." I'm not sure how complex my libaries are and whether they are sufficient or not. 

I think I will move forward with Michael's pipeline for now since he adapted it from Misha Matz's lab, which developed the 3' TagSeq protocol for corals at UT Austin. 
