---
layout: post
title: QC scripts for Pegasus Stress-Hardening RNA-Seq Experiment
date: '2023-05-24'
categories: coding
tags: [coding]
---

These scripts are what I want to use for QC once I get access to Pegasus and can move all my sequences onto the project space. 

#!/bin/bash
#/ this second line in the job script is the relative path, then the next line is the absolute path
```{bash}

#!/bin/bash
#~/scripts/fastqc/fastqc_stresshardening2022.job 
#/scratch/projects/and_transcriptomics/{foldername}/scripts/fastqc/fastqc_stresshardening2022.job
#purpose: quality checking of raw RNAseq reads using FASTQC on Pegasus compute node

#BSUB -J fastqc
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o fastqc_stresshardening2022.out
#BSUB -e fastqc_stresshardening2022.err
#BSUB -n 8
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics/{foldername}" 

cd ${and}
for SAMP in *.fastq.gz

do

module load java/1.8.0_60
module load fastqc/0.10.1 \
${and}/$SAMP \
--outdir ${and}/fastqc_results/
done
```

I've tried running this job several times now (with the proper file folder names) and i keep getting errors related to the module files and the --outdir being unrecognized (probably because fastqc isn't being loaded properly). But I don't know how else to load modules that are already installed on Pegasus, what i wrote in line 36 is what has worked in the past.

Ok i got the most recent version of fastqc_AcerCCC.sh to work! Or at least it's running on Pegasus right now and hasn't immediately failed.

```{bash}
#!/bin/bash
#~/scripts/fastqc_AcerCCC.sh
#/scratch/projects/and_transcriptomics/Allyson_CCC/scripts/fastqc_AcerCCC.sh
#purpose: quality checking of raw RNAseq reads using FASTQC on Pegasus compute node

#BSUB -J AcerCCC_fastqc
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o fastqc_AcerCCC.out
#BSUB -e fastqc_AcerCCC.err
#BSUB -n 8
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics/Allyson_CCC"

cd ${and}
fastqc *.fastq.gz
--outdir ${and}/fastqc/
```

