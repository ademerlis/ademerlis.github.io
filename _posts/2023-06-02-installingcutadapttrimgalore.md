---
layout: post
title: installing cutadapt and trim_galore on Pegasus
date: '2023-06-02'
categories: coding
tags: [coding, CCC_ch4]
---

trim_galore is a wrapper that requires cutadapt and fastQC (https://github.com/FelixKrueger/TrimGalore)

first install cutadapt on pegasus

```{bash}
module load py-pip/20.2
pip install cutadapt
```

That worked.

Now install trim_galore

```{bash}
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
```

I installed Trim_galore into my programs folder in the and_transcriptomics project space.

Here is the code I am currently running on Pegasus: (got from Natalia Andrade's code: https://github.com/China2302/SCTLD_RRC/blob/main/hpc/trimming.sh)

```{bash}
#!/bin/bash
#BSUB -J trim_CCC
#BSUB -q bigmem
#BSUB -P and_transcriptomics
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -o trim_CCC%J.out
#BSUB -e trim_CCC%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

module load fastqc/0.10.1

for sample in ${and}/Allyson_CCC/fastq_files/*.gz ;

do \

${and}/programs/TrimGalore-0.6.10/trim_galore ${sample}
--fastqc \
--fastqc_args "--outdir ${and}/Allyson_CCC/trimmed/" \
--illumina \
--cores 4 \
--three_prime_clip_R1 12 \
--nextseq 30 \
--length 20 \
-o ${and}/Allyson_CCC/trimmed/ ; \

done

multiqc ${and}/Allyson_CCC/trimmed/ \
--outdir ${and}/Allyson_CCC/trimmed/
```

When I submit this job, though, it says that "multicore support is not enabled. proceeding with single-core trimming." Which is annoying because i submitted it to the "bigmem" queue. 

Ok, I found on Pegasus docs that your project has to be given permission to submit to the bigmem queue. I can, however, submit to the parallel queue if I use a specific flag. Parallel is for 16 or more cores, which is what I specified in my job anyways.

The following notes came up when I ran trim_galore:
PLEASE NOTE: Using multi-cores for trimming with 'gzip' only has only very limited effect! (see here: https://github.com/FelixKrueger/TrimGalore/issues/16#issuecomment-458557103)
To increase performance, please install 'pigz' and run again

Proceeding with 'gzip' for decompression
To decrease CPU usage of decompression, please install 'igzip' and run again

I am now trying to install pigz, which I can't use "sudo install" because I am not an administrator on Pegasus. I then tried loading the anaconda environment, but this didn't work either (I don't have permissions to write files to the anaconda shared folder):
```{bash}
source /share/apps/anaconda/anaconda3_build/bin/activate
conda install -c conda-forge pigz
```
EnvironmentNotWritableError: The current user does not have write permissions to the target environment.
  environment location: /share/apps/anaconda/anaconda3_build

Ok i was able to install pigz by downloading the .tar.gz file from here (https://zlib.net/pigz/) and then moving it to Pegasus using scp, and then running "tar xvzf pigz-2.7.tar.gz" in the programs folder.

Ok now I'm going to load pigz in the job script and then also specific the parallel queue instead of bigmem and see if that works faster.

Loading pigz and fastQC within the for loop does NOT work.

```{bash}
and="/scratch/projects/and_transcriptomics"

for sample in ${and}/Allyson_CCC/fastq_files/*.gz ;

do \
module load fastqc/0.10.1 \
${and}/programs/pigz-2.7 \
${and}/programs/TrimGalore-0.6.10/trim_galore ${sample} \
```
I get errors for everything, nothing worked.

I found in the notes that running multi-cores (i.e. using pigz) to unzip fastq files using cutadapt requires Python3, and when I run my script, it says that it cannot detect the version of Python available so I'm not sure if it's able to run pigz even though it's installed now. I tried to run "module load python/3.8.7" within the for loop of the script, and it did not like that. I don't think you can load packages that way. I tried to run "module load python/3.8.7" outside the for loop too, but that didn't do anything.

I also can't get my job to get submitted to the parallel queue even though i used the flags they said to specify on the Pegasus website (#BSUB -q parallel
#BSUB -R "span[ptile=16]"). 

I might just give up and let it run as-is for now, even though it's really slow. I'll see if it times out or if it makes it through the weekend.
