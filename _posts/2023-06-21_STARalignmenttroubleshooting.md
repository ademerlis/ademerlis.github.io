---
layout: post
title: STAR alignment troubleshooting
date: '2023-06-21'
categories: coding
tags: [coding, CCC_ch4, temperaturevariability2023]
---

I'm looking at Young et al. 2020 methodology because they used Iliana Baums' Apal genome, which is also available on Galaxy, so maybe the code required will be similar. 

First, it looks like they installed a local version of the STAR program instead of using the one available on Pegasus. I can't find the version of STAR that was used for their analysis, I guess whatever was available in 2020? 

Ok, a couple more things to note. I am currently looking at their STAR index code for the Apal genome, and I noticed that they 1) removed tRNA from the gff3 file and 2) used the gff3 file in place of the GTF file. (see below)

```{bash}
#!/bin/bash
#BSUB -J STAR_index
#BSUB -q bigmem
#BSUB -P transcriptomics
#BSUB -n 8
#BSUB -R "rusage[mem=4000]"
#BSUB -e /scratch/projects/transcriptomics/ben_young/apalm_v2/scripts/error_and_outputs/STARindex.e%J
#BSUB -o /scratch/projects/transcriptomics/ben_young/apalm_v2/scripts/error_and_outputs/STARindex.o%J

########### MAKE SURE TO GET RID OF tRNA FROM THE GFF FILE
## removing the tRNA lines from the gff3 file
awk '$3 != "tRNA" {print $0}' < /scratch/projects/transcriptomics/ben_young/apalm_v2/apal_genome/v2/Apalm_assembly_v2.0_180910.gff3 > /scratch/projects/transcriptomics/ben_young/apalm_v2/apal_genome/v2/Apalm_assembly_v2.0_180910_no_tRNA.gff3

/nethome/bdy8/Ben_Xaymara_GE_project/programs/STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir /scratch/projects/transcriptomics/ben_young/apalm_v2/star_index/ \
--genomeFastaFiles /scratch/projects/transcriptomics/ben_young/apalm_v2/apal_genome/v2/Apalm_assembly_v2.0_180910.fasta \
--sjdbGTFfile /scratch/projects/transcriptomics/ben_young/apalm_v2/apal_genome/v2/Apalm_assembly_v2.0_180910_no_tRNA.gff3 \
--sjdbOverhang 100 \
--sjdbGTFtagExonParentTranscript Parent
```

I googled whether it is necessary to remove tRNA, because I've never seen that before, I don't think it is necessary.

<img width="688" alt="Screen Shot 2023-06-21 at 11 06 37 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/fafeb23d-153d-4a0a-aabe-eca245ee6ffd">

source: https://biohpc.cornell.edu/doc/RNA-Seq-2018-Lecture1.pdf 

From the same source, I also found this which was interesting (a new tool to convert gff3 to gtf):

<img width="743" alt="Screen Shot 2023-06-21 at 11 22 08 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/6a49733a-ff68-4384-9245-0428b21bdef2">

I then found this website which explains the difference between GFF and GTF and how you need to "adjust them for the intended reading tool." https://goenomics.com/glossary.html
<img width="852" alt="Screen Shot 2023-06-21 at 11 24 09 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/c63a0fd0-fd65-4574-a68e-d3fac3374dc6">

so, first I need to check if I indexed the Acer genome in the same way. That code worked for me. 

