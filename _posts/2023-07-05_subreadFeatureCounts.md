---
layout: post
title: Subread and FeatureCounts
date: '2023-07-05'
categories: coding
tags: [coding, CCC_ch4]
---

I am trying to see if I can still get gene counts on the Acer CCC samples without using StringTie, since that has been difficult for me to get to work. I adapted this featureCounts code from my wound-healing project and submitted the job to pegasus:

```{bash}
#!/bin/bash
#BSUB -J featurecounts_trimmed
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o featurecounts%J.out
#BSUB -e featurecounts%J.err
#BSUB -n 8

and="/scratch/projects/and_transcriptomics"

module load subread/1.6.2

featureCounts -t gene \
-g ID \
-a ${and}/genomes/Acer/Acerv_assembly_v1.0.gff3 \
-o ${and}/Allyson_CCC/subread_counts/AcerCCC.counts \
${and}/Allyson_CCC/aligned/*Aligned.sortedByCoord.out.bam
```

The program ran, but I got very low assigned read rates (also known as alignment rate)of less than 25%. ChatGPT said that alignment rates of 70-90% are good, so something is going on here. 

I am going to try to redo this step using the updated annotation files and see if I get higher alignment rates.

```{bash}
#!/bin/bash
#BSUB -J featurecounts_trimmed_updatedannotations
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o featurecounts_updatedannotations%J.out
#BSUB -e featurecounts_updatedannotations%J.err
#BSUB -n 8

and="/scratch/projects/and_transcriptomics"

module load subread/1.6.2

featureCounts -t gene \
-g ID \
-a ${and}/genomes/Acer/Acerv.GFFannotations.fixed_transcript.gff3 \
-o ${and}/Allyson_CCC/subread_counts/AcerCCC_fixedannotations.counts \
${and}/Allyson_CCC/aligned_updatedannotations/*Aligned.sortedByCoord.out.bam
```

It looks like I'm getting different alignment rates, but still low (less than 25%).

So I loaded the subread module that is on Pegasus to run this, and the version 1.6.2 was from 2019. ChatGPT said that the version of Subread can also significantly affect alignment rates, so maybe I should download the most up-to-date version to my scratch space and try running that instead.

```{bash}
Allysons-MacBook-Pro-2:Downloads allysondemerlis$ scp subread-2.0.6-Linux-x86_64.tar.gz and128@pegasus.ccs.miami.edu:/scratch/projects/and_transcriptomics/programs

[and128@login4 programs]$ tar -xzf subread-2.0.6-Linux-x86_64.tar.gz
cd subread-2.0.6-Linux-x86_64/bin/
nano ~/.bash_profile
#in bash profile add: export PATH=$PATH:/scratch/projects/and_transcriptomics/programs/subread-2.0.6-Linux-x86_64
```

Now I should be able to run the scripts with this program. 

We'll update the original script first:
```{bash}
#!/bin/bash
#BSUB -J featurecounts_trimmed
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o featurecounts%J.out
#BSUB -e featurecounts%J.err
#BSUB -n 8

and="/scratch/projects/and_transcriptomics"

/scratch/projects/and_transcriptomics/programs/subread-2.0.6-Linux-x86_64/bin/featureCounts -t gene \
-g ID \
-a ${and}/genomes/Acer/Acerv_assembly_v1.0.gff3 \
-o ${and}/Allyson_CCC/subread_counts/AcerCCC.counts \
${and}/Allyson_CCC/aligned/*Aligned.sortedByCoord.out.bam
```

It looks like I am getting the same percentages. Ugh ok whatever. I'll just run it on the updated annotations files too just so everything is updated.


```{bash}
#!/bin/bash
#BSUB -J featurecounts_trimmed_updatedannotations
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o featurecounts_updatedannotations%J.out
#BSUB -e featurecounts_updatedannotations%J.err
#BSUB -n 8

and="/scratch/projects/and_transcriptomics"

/scratch/projects/and_transcriptomics/programs/subread-2.0.6-Linux-x86_64/bin/featureCounts -t gene \
-g ID \
-a ${and}/genomes/Acer/Acerv.GFFannotations.fixed_transcript.gff3 \
-o ${and}/Allyson_CCC/subread_counts/AcerCCC_fixedannotations.counts \
${and}/Allyson_CCC/aligned_updatedannotations/*Aligned.sortedByCoord.out.bam
```

NOTE: you did not change the names of the output files so you overwrote them. whoops. but the only difference is the version of subread and it didn't look like it changed percent alignment very much (still around 25%).
