---
layout: post
title: md5 check sums script
date: '2023-07-19'
categories: coding
tags: [coding, temperaturevariability2023, CCC_ch4]
---

So I followed the first script from [Sam](https://github.com/SamGurr/SamGurr.github.io/blob/master/_posts/2021-01-07-Geoduck-TagSeq-Pipeline.md) and [Ariana](https://github.com/AHuffmyer/EarlyLifeHistory_Energetics/blob/master/Mcap2020/Scripts/TagSeq/Genome_V3/TagSeq_BioInf_genomeV3.md) to do the "md5 check sums" thing. I guess the purpose of it is to check the integrity of the files after downloading them from Basespace. Idk why it's necessary but everyone in Hollie's lab seems to do it, so it must be important! However, I submitted this job on Pegasus yesterday for the Ch2_tempvariability2023 raw fastq.gz files and the job had still not yet started today:

```{bash}
#!/bin/bash
#BSUB -J transfer_checks
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o transfer_checks.out
#BSUB -e transfer_checks.err
#BSUB -n 8
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023"

# generate md5
md5sum ${and}/fastq_rawreads/*.gz > ${and}/tempvariability.md5

# count number of sequences per file

zcat ${and}/fastq_rawreads/*.fastq.gz | echo $((`wc -l`/4)) > ${and}/fastq_rawreads/rawread.counts.txt
```

I'm not sure what's wrong with it and why it hasn't started. Does it require a lot of memory or something? In Sam code, this was what they specified in the job:
```{bash}
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH -D /data/putnamlab/KITT/hputnam/20201217_Geoduck_TagSeq/
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3
```

In Ariana's code, it looks like she didn't even submit a script, just ran it in the terminal? 

<img width="1040" alt="Screen Shot 2023-07-19 at 10 16 43 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/68119471-0dad-4f84-80d1-ad0f8d4c5d44">

In [Zoe's](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/ZD_Heron-Pdam-gene-expression.md), she doesn't even have the md5 check sums code, and her's is the most recent pipeline, so I'm just going to skip this step.


