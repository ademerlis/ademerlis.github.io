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

None of the above workflows specify if there is a min/max amount of sequence files to use. Dr. Jill Ashey used 12 samples (6 forward, 6 reverse). I emailed Jill and she recommended: 
> "I think that for Trinity, the number of samples isn't as important as the number of reads. Typically, RNAseq samples have about 20-50M reads per sample. If we are going off the Trinity recommendation that 50M is saturated, we would only need a few samples for the transcriptome. When creating a de novo transcriptome for a species that doesn't have any genomic resources, I personally think that it's better to include more samples across a range of conditions (n=3-5 per condition) so that more transcripts can be captured across the board and a proper meta-transcriptome can be assembled."

From [Grabherr et al. 2011 (the first Trinity publication)](https://www.nature.com/articles/nbt.1883#MOESM13), they state that 50 M pairs of reads is enough for Trinity to fully reconstruct 86% of annotated transcripts. The authors say that after 50 M paired reads, it is "saturated" or doesn't require any more data. 

[Haas et al. 2014 (the newer Trinity publication)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3875132/) doesn't talk about this. 


## *de novo* transcriptome assembly publications 

1. [Raghavan et al. 2022](https://doi.org/10.1093/bib/bbab563) -
   - this paper has good conceptual figures
   - recommends in-silico normalization
   - Trinity is De Bruijn graph-based 
   - post-assembly QC is important
      - N50 value is a common sequence length statistic
      - however, N50 value should be used with caution - the ExN50 is a better modification for transcriptome assemblies because the recovery has been for short full-length sequences, rather than few, very long contigs
      - ExN50 = expression-weighted sum of corresponding isoform lengths
      - this metric is only implemented for Trinity 
      - the fraction of all reads that map back to the assembly 
      - proportion of reads that map to multiple sequences should be low
      - BUSCO (Benchmarking Universal Single-Copy Orthologs) is another good tool to determine if a large majority of orthologs are found to the assembled transcriptome. Looking for "universal" genes. General rule is a good BUSCO completeness score is >80%
      - Tools for QC include TransRate, DETONATE, and rnaQUAST


## Coral papers that construct *de novo* transcriptomes

1. [Alderdice et al. 2022](https://www.nature.com/articles/s41598-022-22604-3) -
   - Paired end reads (150bp) from *Acropora*
   - SOAPdenovo-Trans for 23-kmer length was used for de novo transcriptome assembly from **all samples**.
   - Number of samples: 32
   - Extraction protocol: Qiagen RNeasy mini kit 
   - Sequencing: NovaSeq 6000

2. [Lock et al. 2022](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.9345) - 
   - Paired end reads (150bp) from *Porites lobata*
   - All reads combined and normalized using in-silico normalization from Trinity
   - All samples/reads used for de novo transcriptome assembly with Trinity
   - Number of samples: 60
   - Extraction protocol: NEBNext Ultra RNA Library Prep Kit (New England Biolabs)
   - Sequencing: NextSeq 500

3. [Poquita-Du et al. 2019](https://doi.org/10.3389/fmars.2019.00121) -
   - *Pocillopora acuta*
   - All sequence reads from the same genotype (so 12 samples?) were combined and assembled de novo using Trinity -- each genotype got their own individual reference assembly
   - Number of samples: 36
   - Extraction protocol: TRIzol following Barshis et al. 2013
   - Sequencing: Illumina-based cDNA library prep + Illumina HiSeq 2500

4. [Poquita-Du et al. 2020](https://doi.org/10.1007/s00338-020-01902-0) -
   - they used the same sequencing data as in #3, but this paper focuses on the symbiont reads (which they obtained from the same de novo assembly)

5. [Anderson et al. 2016](http://dx.doi.org/10.7717/peerj.1616) -
   - *Orbicella faveolata*
   - TRIzol extractions
   - Illumina GAIIx sequencing
   - Trinity used for de novo assembly
   - 6 samples were used, with 3 of them being from diseased coral and 1 being from a bleached coral

6. [Veglia et al. 2018](https://doi.org/10.1016/j.margen.2018.08.003 - 
   - high quality reference genome of *Agaricia lamarcki*
   - TRIzol extraction
   - Total RNA sequenced on Illumina Highseq4000 for poly-A tail seelction and 150bp paired-end reads
   - They used a de novo transcriptome assembly pipeline (http://github.com/NCGAS/de-novo-transcriptome-assembly-pipeline) that tests all the different assembly softwares (Trinity, SOAPdenovo-Trans, Velvet, Oases, Trans-ABySS) - combines all those assemblies together. Then, a "clean consensus" transcriptome was generated using EvidentialGene.
   - Looks like it is just 1 tissue sample




   



