---
layout: post
title: Pcli Salmon transcriptome
date: '2023-07-26'
categories: coding
tags: [coding, Ch2_tempvariability]
---

So I installed Salmon locally onto my scratch space by downloading the Linux binary package:

```{bash}
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz
tar -xzvf salmon-1.5.2_linux_x86_64.tar.gz
cd salmon-1.5.2_linux_x86_64/bin/
#salmon is green so that means it is already compiled
#add to PATH
nano ~/.bash_profile
source ~/.bash_profile
#check it works:
salmon-1.5.2_linux_x86_64/bin/salmon --verson
```

This salmon code worked:
```{bash}
#!/bin/bash
#BSUB -J Pcli_transcriptome_salmon_index
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -o Pcli_transcriptome_salmon_index%J.out
#BSUB -e Pcli_transcriptome_salmon_index%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

${and}/programs/salmon-1.5.2_linux_x86_64/bin/salmon index -t ${and}/genomes/Pcli/clean_Pcli_transcriptome_final.fasta -i ${and}/genomes/Pcli/Pcli_transcriptome_index
```

Output:
![Screen Shot 2023-07-26 at 9 17 51 AM](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/1acc45ef-e328-483d-a193-6ba5535cc124)

Following [this Salmon tutorial](https://combine-lab.github.io/salmon/getting_started/):

```{bash}
#!/bin/bash
#BSUB -J Pcli_transcriptome_salmon_quant
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -o Pcli_transcriptome_salmon_quant%J.out
#BSUB -e Pcli_transcriptome_salmon_quant%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/AS_pipeline/3_trimmed_fastq_files/Pcli_fastq_files"

data=($(ls *.gz))

for sample in ${data[@]} ;

do \
${and}/programs/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${and}/genomes/Pcli/Pcli_transcriptome_index -l A -1 ${sample} -p 8 --validateMappings -o ${and}/Ch2_temperaturevariability2023/AS_pipeline/4_Pcli_specific/salmon_quant_files ; \
done
```
