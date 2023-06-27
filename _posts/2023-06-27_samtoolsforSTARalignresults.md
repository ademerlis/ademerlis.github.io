---
layout: post
title: samtools for STAR alignment results
date: '2023-06-27'
categories: coding
tags: [coding, CCC_ch4]
---

So I went down the path of using samtools after STAR alignment because it seems that it is used as a quality-checking tool to see if alignment worked. I thought multiQC did that but idk what the difference is. So I tried first running Danielle's code (https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-04-14-Molecular-Underpinnings-RNAseq-Workflow.md) but then realized her sequences are paired-end reads and that's why the code wasn't working. 

Then I tried running this code and it didn't work either (it was like stuck in the job run mode but wasn't outputting anything so I had to kill the job):

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
samtools flagstat ${i}
done
```

But then I found this biostar thread online that had different samtools commands (https://www.biostars.org/p/138116/) and I ran this and it "worked" but didn't print the sample info with it:
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
samtools view -F 0x904 -c ${i}
done
```

The output file printed this:
```{bash}
451173
14339859
14922510
12474659
11287053
9161911
3417981
2868510
7828366
10810654
10619206
13883067
14504107
12929727
6133175
391307
6845324
```

Based on the biostar thread, "samtools view -F 0x904 -c file.bam" should output the number of mapped alignments, excluding multimappers for a single-end dataset.

But when I look at the multiqc results of the STAR_align_trimmed samples, I see that the numbers from above correspond to the number of uniquely mapped + mapped to multiple loci (image for example below):

![Screen Shot 2023-06-27 at 9 09 16 AM](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/5f12614b-155c-4d6e-8854-4811dfe03340)

So idk what's better, samtools of multiqc. And does it really matter? How do I verify the quality of the alignment...

Also here's a meme to summarize this week:

![FC72BDD9-EFF5-44D3-A5B8-3F269F3F283B](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/81988755-537e-4998-a9cf-a185ff7e2a66)
