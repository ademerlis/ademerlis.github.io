---
layout: post
title: Starting over
date: '2023-07-18'
categories: coding
tags: [coding, temperaturevariability2023, CCC_ch4]
---

I think I need to start from the beginning and go through one pipeline that has been proven for that person (or lab group) to work time and time again, rather than try to frankenstein pieces of people's codes together to get something to work (which is what I have tried thus far and haven't been successful). It also helps to find a person with well-annotated code. Thankfully, Hollie Putnam's lab has tried-and-true methods that go back to pipelines of [Dr. Sam Barr](https://github.com/SamGurr/SamGurr.github.io/blob/master/_posts/2021-01-07-Geoduck-TagSeq-Pipeline.md) and [Dr. Ariana Huffmyer](https://github.com/AHuffmyer/EarlyLifeHistory_Energetics/blob/master/Mcap2020/Scripts/TagSeq/Genome_V3/TagSeq_BioInf_genomeV3.md). 

The reason I am doing this stems from discussions I had both with Jill Ashey (who pointed me towards the most recent pipeline from Hollie's lab, [Zoe Dellaert](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/ZD_Heron-Pdam-gene-expression.md)), and [Kevin Wong](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/scripts/TagSeq/TagSeq_Analysis_HPC.md) (who made a good point to me that a TagSeq pipeline and an RNAseq pipeline are very different things, and that I can't frankenstein codes from people doing either method). 

Thus, I will follow the TagSeq pipelines of Dr. Barr and Dr. Huffmyer, and follow any updates that Zoe made in her code since she is the most recent to do it (I think). 

First, let me review what I have attempted to do so far to give a summary for future me:

I was trying to analyze the Ch2_tempvariability Acer samples and the Ch4_AcerCCC samples the same way (but Ch 2 was TagSeq from UT Austin, and Ch 4 was QuantSeq (which equals RNAseq??) with UM Genomics).

Pipeline included:
1. TrimGalore for adapter (from [Natalia's code](https://github.com/China2302/SCTLD_RRC/blob/main/hpc/trimming.sh)
2. TrimGalore for polyA tail (from Natalia's code as well)
3. Create STAR index (following a combo of [Natalia's code](https://github.com/China2302/SCTLD_RRC/blob/main/hpc/STAR_index.sh) and [Jill's code](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/RNASeq_pipeline_FL.md)
4. STAR alignment (following a combo of [Natalia's code](https://github.com/China2302/SCTLD_RRC/blob/main/hpc/star_align_trimmed.sh) and [Jill's code](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/RNASeq_pipeline_FL.md)
5. Stringtie + merge GTF + gffcompare + re-assemble with Stringtie (following [Jill's code](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/RNASeq_pipeline_FL.md)
6. for the Acer CCC I also tried to use featurecounts following my code from the [Pdam wound healing 2019](https://github.com/ademerlis/woundhealingPdam/blob/main/bioinformatics/6_featureCounts.txt) and tried to quantify directly from STAR following [Natalia's code](https://github.com/China2302/SCTLD_RRC/blob/main/05_read_counts.Rmd)

Important to note that I don't think Natalia or Jill used TagSeq... so that could be why mine isn't working (at least for the Ch2_tempvariability one. Who knows about the C4_AcerCCC). 

What I should try (following Sam, Zoe, and Ariana's pipelines):
1. md5 to check file integrity
3. Trimming with fastp
4. alignment with HISAT2
5. assembly with Stringtie2
6. generate counts matrix with prepDE.py from Stringtie2
