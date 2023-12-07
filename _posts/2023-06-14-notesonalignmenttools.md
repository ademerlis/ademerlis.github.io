---
layout: post
title: notes on alignment tools
date: '2023-06-14'
categories: coding
tags: [coding, CCC_ch4, temperaturevariability2023]
---

 In the past, for the Pdamicornis wound healing paper, I used STAR to align the 3' RNAseq reads to the genome, not the transcriptome. 
 
 The notes I'm relying on for this round of analysis are from Dr. Natalia Andrade (https://github.com/China2302/transcriptomics-workflow) and Dr. Michael Studivan (https://github.com/mstudiva/tag-based_RNAseq/blob/master/). The reason for this is because I sequenced the Acer CCC samples with Natalia's SCTLD Mote and Nova samples, while the temperature variability 2023 Acer and Pcli samples were sequenced following Michael's protocols. 
 
 So far, because Michael's stuff is all in python, I have had more success with Natalia's pipeline. However, the QC and trimming steps are not super specific, and they both used cutadapt to trim Illumina adapters. Michael's code had more specifics for trimming (i.e. specific base pair sequences) that are potentially more likely to be found in TagSeq generated reads. Natalia specifically trimmed the polyA tails because that is a problem in 3' RNA seq. 
 
 But for alignment, there are many different possible tools to use, and there is also the question of whether to align to a genome or a transcriptome. 
 
 The UM CCS student mentors made a pros and cons list for aligning to a genome versus a transcriptome (https://github.com/ccsstudentmentors/tutorials/tree/master/RNA-Seq/Quantifying-RNA-Expression). I think it also depends on what is available for the species you're analyzing. 
 
 Transcriptome -- then use Bowtie.
 
 Genome -- use STAR.
 
 After talking to Dr. Kevin Wong, he said to just pick whichever gives the best alignment rate. 
 
 In terms of which tool to use, STAR seems to be the superior aligner:
 1. https://www.biostars.org/p/353946/
 2. <img width="737" alt="Screen Shot 2023-06-15 at 8 17 33 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/785896b7-c784-48bd-9d2b-995030d55f94">

Also important to note for RNA-Seq (with polyA+ and total) is using an aligner that "looks for splicing junctions in addition to genome mapping". STAR is a fast splice read aligner that uses short reads. The key limitation is RAM. (http://homer.ucsd.edu/homer/basicTutorial/mapping.html)

To use STAR (notes come from this tutorial: http://homer.ucsd.edu/homer/basicTutorial/mapping.html):

1. Build a genome index (common step for all aligners). Make a directory for the index, then copy the genome FASTA files into the directory. Make that your current directory and then run a script (or submit a job on an HPC) to build the index. 

```{bash}
STAR  --runMode genomeGenerate --runThreadN <# cpus> --genomeDir <genome output directory> --genomeFastaFiles <input Genome FASTA file>
```

2. Align RNA-Seq reads to the genome with STAR


```{bash}
STAR --genomeDir <Directory with the Genome Index>  --runThreadN <# cpus> --readFilesIn <FASTQ file> --outFileNamePrefix <OutputPrefix>
#note: this is for single-end data
```

STAR will create several output files. The ".aligned.out.sam" file is the most important. Convert it to BAM file to use with other programs.

