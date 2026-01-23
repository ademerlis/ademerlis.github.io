---
layout: post
title: Next steps with PSTR RT transcriptomics
date: '2026-01-08'
categories: Analysis
tags: [reciprocal transplant, coding, Pegasus, rRNA, de novo transcriptome assembly]
---

Some of the trimmed sequence files (<2GB) were run on sortmeRNA using --blast, and the others were run more recently (using [the most up to date sortmeRNA code](https://github.com/ademerlis/reciprocaltransplant/blob/main/gene_expression/6_sortmerna_redo.sh)), which did not include --blast.

This resulted in differing files (no singletons) -- asking chatGPT, they said that sortmeRNA doesn't change the filtering parameters for rRNA vs non_rRNA, it just doesn't pull out paired vs singletons.

So in essence, I can use the "old" files and the "new" files together to assemble a transcriptome using trinity.

creating Trinity environment in conda: 
`conda create -n trinity -c conda-forge -c bioconda trinity`

moved all the _non_rRNA and _rRNA files up, and then identified the _out_paired from the first round and renamed those to non_rRNA and moved those too

now identify samples to use for transcriptome assembly (i selected a subset across genotypes and time points, and then looked for ones that had fwd and rev sortmerna files > 1 GB but less than 2 GB).

look for samples from sortmerna with < 1 GB -- rerun those (aka the one at 17 k)

samples to rerun (all at 17kb):
- Pstr_Jul2023_006
- Pstr_Nov2023_180
- Pstr_Nov2023_181
- Pstr_Nov2023_190
- Pstr_Nov2023_191
- Pstr_Nov2023_195
- Pstr_Nov2023_196

Jan 9, 2026: I reran sortmeRNA on these and it worked, so now they are all ready to go.

This is what ended up working to install and activate trinity on pegasus (although some warnings and errors came up that may cause problems later):
```{bash}
conda update trinity
#to get mamba to work on pegasus:
module load mambaforge/1.5.8
source /share/apps/mambaforge/install/etc/profile.d/conda.sh
source /share/apps/mambaforge/install/etc/profile.d/mamba.sh
mamba install trinity
mamba update trinity
conda activate /nethome/and128/anaconda3/envs/trinity
which Trinity
conda env create --name trinity --file /nethome/and128/anaconda3/envs/trinity/environment.yml
conda env export > environment.yml
conda env create --name trinity --file environment.yml
```

I ran this, and I keep running into errors:
```{bash}
command="Trinity --seqType fq --left pstr_fwd.fq.gz --right pstr_rev.fq.gz --CPU 10 --max_memory 100G --min_kmer_cov 2 --output /scratch/projects/and_transcriptomics/reciprocaltransplant/raw_seq_files/trinity/output"
bsub -P and_transcriptomics -q bigmem -n 10 -R "rusage[mem=10000]" -W 120:00 -J trin_pstr -e trin_pstr.err -o trin_pstr.out eval ${command}
```

```{bash}
#!/bin/bash
#BSUB -P and_transcriptomics
#BSUB -q bigmem
#BSUB -n 10
#BSUB -R "rusage[mem=10000]"
#BSUB -W 120:00
#BSUB -J trin_pstr
#BSUB -e trin_pstr.err
#BSUB -o trin_pstr.out

# from pegasus docs
module load miniforge3/24.3.0-0

# Activate Trinity environment
conda activate /nethome/and128/anaconda3/envs/trinity

# Run Trinity
Trinity --seqType fq \
  --left /scratch/projects/and_transcriptomics/reciprocaltransplant/sortmerna/trinity/pstr_fwd.fq.gz \
  --right /scratch/projects/and_transcriptomics/reciprocaltransplant/sortmerna/trinity/pstr_rev.fq.gz \
  --CPU 10 \
  --max_memory 100G \
  --min_kmer_cov 2 \
  --output /scratch/projects/and_transcriptomics/reciprocaltransplant/raw_seq_files/trinity/output
```

The most recent errors are:
error while loading shared libraries: libncurses.so.5: cannot open shared object file: No such file or directory
This Perl not built to support threads
Compilation failed in require at /nethome/and128/anaconda3/envs/trinity/bin/Trinity line 5.
BEGIN failed--compilation aborted at /nethome/and128/anaconda3/envs/trinity/bin/Trinity line 5.

So i am trying this:
`conda install -c bioconda -c conda-forge trinity samtools>=1.10 perl-threaded`

It still is having issues. I wonder if it is because I have another conda environment with BioPerl. 


**Jan 23 2026**

I was able to run Trinity in the trinity conda environment, and got the job submitted using this directly in the command line:
```{bash}
command="Trinity --seqType fq --left /scratch/projects/and_transcriptomics/reciprocaltransplant/raw_seq_files/trinity/pstr_fwd.fq.gz --right /scratch/projects/and_transcriptomics/reciprocaltransplant/raw_seq_files/trinity/pstr_rev.fq.gz --CPU 10 --max_memory 100G --min_kmer_cov 2 --output /scratch/projects/and_transcriptomics/reciprocaltransplant/raw_seq_files/trinity/trinity"
bsub -P and_transcriptomics -q bigmem -n 10 -R "rusage[mem=10000]" -W 120:00 -J trin_pstr -e trin_pstr.err -o /scratch/projects/and_transcriptomics/reciprocaltransplant/raw_seq_files/trinity/output eval ${command}
```

However, it quit after a bit and gave me 2 output files, which basically had the error message:
`Argument "error" isn't numeric in numeric ne (!=) at /nethome/and128/anaconda3/envs/trinity/bin/Trinity line 3812.
Error, need samtools installed that is at least as new as version 1.3 at /nethome/and128/anaconda3/envs/trinity/bin/Trinity line 3813.`

Now i'm trying to install samtools again and I'm running into the same issue as before.

