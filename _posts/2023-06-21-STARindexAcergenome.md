---
layout: post
title: STAR index Acer genome
date: '2023-06-21'
categories: coding
tags: [coding, CCC_ch4, temperaturevariability2023]
---

The first step when aligning reads to the genome when using STAR is to first create an index based on the most up-to-date genome assembly. I downloaded a 2019 draft assembly from Baums et al. (https://usegalaxy.org/published/history?id=1f8678b27ae56467 and https://usegalaxy.org/published/page?id=2f8d21c73f8501e2 for file descriptions based on Apal files). 

When I indexed the genome, I used this fasta file: "Acerv_assembly_v1.0_171209.fasta" which contains everything. There are also .fa files available which had just mRNA vs protein. there is also a "masked.fa" file, which I don't know what that means (on the Galaxy page it says "repeat-masked scaffolds"). 

I think I can just use the fasta file I already used. I made a couple updates to the index script, so this is the most up-to-date version:
```{bash}
#!/bin/bash
#BSUB -J Acer_star_index
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -o /scratch/projects/and_transcriptomics/Genomes/Acer/star_index%J.out
#BSUB -e /scratch/projects/and_transcriptomics/Genomes/Acer/star_index%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

module load star/2.3.0e

STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeDir ${and}/genomes/Acer/Acer_STAR_index \
--genomeFastaFiles ${and}/genomes/Acer/Acerv_assembly_v1.0_171209.fasta \
--sjdbGTFfile ${and}/genomes/Acer/Acerv_assembly_v1.0.gff3 \
--sjdbOverhang 100 \
--sjdbGTFtagExonParentTranscript Parent
```

I combined my original code with updates based on Young et al 2020 and the STAR manual (STAR manual says for the --sjdbOverhang flag, a default of 100 should work for most cases). 

Update: This is the most up-to-date index code that worked:

```{bash}
#!/bin/bash
#BSUB -J Acer_star_index
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -o /scratch/projects/and_transcriptomics/Genomes/Acer/star_index%J.out
#BSUB -e /scratch/projects/and_transcriptomics/Genomes/Acer/star_index%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

/scratch/projects/and_transcriptomics/programs/STAR-2.7.10b/bin/Linux_x86_64_static/STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeDir ${and}/genomes/Acer/Acer_STAR_index \
--genomeFastaFiles ${and}/genomes/Acer/Acerv_assembly_v1.0_171209.fasta \
--sjdbGTFfile ${and}/genomes/Acer/Acerv_assembly_v1.0.gff3 \
--sjdbOverhang 100 \
--sjdbGTFtagExonParentTranscript Parent
--genomeSAindexNbases 13
```


