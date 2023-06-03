---
layout: post
title: writing a script to create multiple jobs at once for trimming
date: '2023-06-03'
categories: coding
tags: [coding, CCC_ch4]
---

This is the code I have right now for the trimming script (note: the parallel flag didn't work, so it ran each file individually one at a time on Pegasus in the general queue, which was very slow):

```{bash}
#!/bin/bash
#BSUB -J trim_CCC
#BSUB -q parallel
#BSUB -R "span[ptile=16]"
#BSUB -P and_transcriptomics
#BSUB -n 16
#BSUB -W 120:00
#BSUB -o trim_CCC%J.out
#BSUB -e trim_CCC%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

for sample in ${and}/Allyson_CCC/fastq_files/*.gz ;

do \
${and}/programs/TrimGalore-0.6.10/trim_galore ${sample} \
--gzip \
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

I want to try something different, based on codes Mike Connelly and Ben Young gave me back in 2019. They wrote a script that then created a job for each sample. So basically manually creating parallel jobs so they can all run at once. This is what it looked like:

```{bash}
#!/bin/bash
#purpose: trimming
#Thanks Mike Connelly and Ben Young!

#BSUB -J trimmomatic_WH
#BSUB -q general
#BSUB -P transcriptomics
#BSUB -o /projects/scratch/transcriptomics/allysondemerlis/scripts/trimmomatic_WH%J.out
#BSUB -e /projects/scratch/transcriptomics/allysondemerlis/scripts/trimmomatic_WH%J.err

#specify variable containing sequence file prefixes and directory paths
and="/projects/scratch/transcriptomics/allysondemerlis"

#cd ${and}/sequences/zippedreads/combinedreads
samples=${and}/sequences/zippedreads/combinedreads/samples.txt

# trimming the files
for SAMP in `cat ${samples}`
do
echo '#!/bin/bash' > /scratch/projects/transcriptomics/allysondemerlis/scripts/${SAMP}_trimming.job
echo '#BSUB -q general' >> /scratch/projects/transcriptomics/allysondemerlis/scripts/${SAMP}_trimming.job
echo '#BSUB -J '$SAMP'' >> /scratch/projects/transcriptomics/allysondemerlis/scripts/${SAMP}_trimming.job
echo '#BSUB -o /scratch/projects/transcriptomics/allysondemerlis/scripts/'$SAMP'_error_trimming.txt' >> /scratch/projects/transcriptomics/allysondemerlis/scripts/${SAMP}_trimming.job
echo '#BSUB -e /scratch/projects/transcriptomics/allysondemerlis/scripts/'$SAMP'_output_trimming.txt' >> /scratch/projects/transcriptomics/allysondemerlis/scripts/${SAMP}_trimming.job

echo 'module load java/1.8.0_60' >> /scratch/projects/transcriptomics/allysondemerlis/scripts/${SAMP}_trimming.job
echo 'module load trimmomatic/0.36' >> /scratch/projects/transcriptomics/allysondemerlis/scripts/${SAMP}_trimming.job

echo '/share/opt/java/jdk1.8.0_60/bin/java -jar /share/apps/trimmomatic/0.36/trimmomatic-0.36.jar \
PE \
-phred33 \
-trimlog /scratch/projects/transcriptomics/allysondemerlis/scripts/'${SAMP}'_trim.log \
'${and}'/sequences/zippedreads/combinedreads/'${SAMP}'_forward.txt.bz2 '${and}'/sequences/zippedreads/combinedreads/'${SAMP}'_reverse.txt.bz2 \
'${SAMP}'_trimmed_pairedfwd.fastq.gz '${SAMP}'_trimmed_unpairedfwd.fastq.gz '${SAMP}'_trimmed_pairedrev.fastq.gz '${SAMP}'_trimmed_unpairedrev.fastq.gz \
ILLUMINACLIP:/share/apps/trimmomatic/0.36/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:80' >> /scratch/projects/transcriptomics/allysondemerlis/scripts/${SAMP}_trimming.job
bsub < /scratch/projects/transcriptomics/allysondemerlis/scripts/${SAMP}_trimming.job
done
```

Note that they used trimmmomatic instead of cutadapt. Need to look into whether there is a real difference in the two. I doubt it if the goal is just to trim the illumina adapters.

