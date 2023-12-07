---
layout: post
title: STAR index and alignment for Pcli
date: '2023-07-25'
categories: coding
tags: [coding, Ch2_tempvariability]
---

Ok so since I have STAR already installed locally on my scratch space, I'm going to try to do the genome index and alignment using STAR for the Pcli transcriptome. ChatGPT says you can just run it the same way as if you were running a genome index.

```{bash}
#!/bin/bash
#BSUB -J Pcli_transcriptome_STAR_index
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -o Pcli_transcriptome_STAR_index%J.out
#BSUB -e Pcli_transcriptome_STAR_index%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

/scratch/projects/and_transcriptomics/programs/STAR-2.7.10b/bin/Linux_x86_64_static/STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeDir ${and}/genomes/Pcli \
--genomeFastaFiles ${and}/genomes/Pcli/clean_Pcli_transcriptome_final.fasta
```

Ok STAR was taking too long so I'm going to try salmon instead.

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

salmon index -t transcriptome.fasta -i transcriptome_index
```

Ok the salmon one isn't running now either. It makes me think that there's something wrong with the Pcli fasta file or something. 
