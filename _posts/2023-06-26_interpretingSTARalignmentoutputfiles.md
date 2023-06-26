---
layout: post
title: interpreting STAR alignment output files
date: '2023-06-26'
categories: coding
tags: [coding, CCC_ch4]
---

Each sample comes with these output files following STAR alignment script:

<img width="368" alt="Screen Shot 2023-06-24 at 11 22 40 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/3f10b651-5d01-4321-8c6c-686f164cc0e5">

But in looking at other people's pipelines (https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/RNASeq_pipeline_FL.md and https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-04-14-Molecular-Underpinnings-RNAseq-Workflow.md) it looks like you need to do more than just MultiQC the aligned reads to make sure if the alignment worked. Otherwise you could have potential errors where things are aligned to the wrong part (see this thread: https://github.com/alexdobin/STAR/issues/613).

I also went down a mini-rabbit hole because I remembered Jill's code had said that the "transcript_id" needed to be added to the Acer gff3 file because STAR wouldn't run without it. But my STAR code ran, and when I look at the gff3 file I downloaded I see this: 

![Screen Shot 2023-06-26 at 2 43 52 PM](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/2408b54d-96ca-4507-b5ee-5a398b646eee)

so it doesn't have the "transcript_id" but it does have "ID" and the STAR alignment code didn't fail. So what happened? 

Ok I fell into another rabbit hole because I started reading this paper on "Alignment and mapping methodology influence transcript abundance estimation" (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8) and "Errors in RNA-Seq quantification affect genes of relevance to human disease" (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0734-x) which basically say that the quantification tool you picks significantly influences number of DEGs. Great. 

But whatever I'm just going to try to install and run "samtools" today, because that seems like an important step to "check the number of mapped reads" (https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2022-12-10-P.verrucosa-RNASeq-Workflow-Host.md).

Is samtools on Pegasus? (module avail)
Yes, there are three versions: samtools/0.1.19, 1.2, and 1.3. 

I'll try 1.3 and run it in my code. (Section of code below adapted from Danielle Becker's code)
I think we want to use the sortedByCoord.out.bam files based on Danielle's code notes (although it isn't explicitly stated)
```{bash}
#!/bin/bash
#BSUB -J samtools_aligned_trimmed
#BSUB -q bigmem
#BSUB -P and_transcriptomics
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -o samtools_aligned_trimmed%J.out
#BSUB -e samtools_aligned_trimmed%J.err
#BSUB -u and128@miami.edu
#BSUB -N

module load samtools/1.3

array1=($(ls /scratch/projects/and_transcriptomics/Allyson_CCC/aligned/*sortedByCoord.out.bam))  
for i in ${array1[@]}; do
samtools view -F 0x4 ${i} | cut -f 1 | sort | uniq | wc -l > mapped_read_counts
done
```

Ok looking at this github page, I can't use the samtools code that Danielle wrote because hers is specified for paired-end reads and I have single-end reads. (https://gist.github.com/davfre/8596159). I need to use different flags than she did.

It looks like in Danielle's code she also experienced issues with the genome assembly file she used, with an issue relating to transcript_id and some other things within the gff3 file itself. This prevented her from performing gene counts with StringTie. She ended up modifying the genome assembly file (see this thread: https://github.com/Putnam-Lab/Lab_Management/issues/11). So this is something that might come up for me when I try to run stringTie.

I am overall very confused with what to do with my STAR alignment outputs (other than run samtools). So maybe I'll just try stringTie next and see what that does. Subread (which has featureCounts) is installed to Pegasus though so maybe I could try that too since it's already there. 

I also found this 2020 paper that did a comparison of all possible RNA-seq bioinformatics pipelines: https://www.nature.com/articles/s41598-020-76881-x

I'm still not 100% sure whether pseudo-alignment or alignment to a genome is better. Also, is alignment to a transcriptome even an option? I guess if I used the Libro et al. 2013 dataset it would be. Otherwise I only have a genome. 
