---
layout: post
title: stringtie on SH samples
date: '2023-07-05'
categories: coding
tags: [coding, temperaturevariability2023]
---

I am now going to try the stringtie scripts on the stress-hardening Acer samples and see if I get better results than the Acer CCC ones. 

I still get the same error: "/projects/lsf_spool/1688568095.27980042.shell: line 21:  : command not found"

this is the script I ran:
```{bash}
#!/bin/bash
#BSUB -J stringtie
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o stringtie%J.out
#BSUB -e stringtie%J.err
#BSUB -u and128@miami.edu
#BSUB -N

module load python/3.8.7

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Allyson_stresshardening_RNAseq/aligned"

data=($(ls *Aligned.sortedByCoord.out.bam))

for i in ${data[@]} ;

do \
/scratch/projects/and_transcriptomics/programs/stringtie-2.2.1/stringtie -G /scratch/projects/and_transc$
done
```

The next step is to merge the stringTie-produced .gtf files

```{bash}
mkdir stringtie_gtf_files
mv *.gtf stringtie_gtf_files/
ls *gtf > acerv_mergelist.txt
cat acerv_mergelist.txt
```

Now I run the merge script:

```{bash}
#!/bin/bash
#BSUB -J stringtie_mergegtf
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o stringtie_mergegtf%J.out
#BSUB -e stringtie_mergegtf%J.err
#BSUB -u and128@miami.edu
#BSUB -N

module load python/3.8.7

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Allyson_stresshardening_RNAseq/aligned/stringtie_gtf_files"

${and}/programs/stringtie-2.2.1/stringtie --merge -p 8 -G ${and}/genomes/Acer/Acerv_assembly_v1.0.gff3 -o stringtie_acerv_merged.gtf acerv_mergelist.txt
```

Next, assess assembly quality:
```{bash}
cd /scratch/projects/and_transcriptomics/Allyson_stresshardening_RNAseq/aligned/stringtie_gtf_files

/scratch/projects/and_transcriptomics/programs/gffcompare-0.12.6.Linux_x86_64/gffcompare -r /scratch/projects/and_transcriptomics/genomes/Acer/Acerv_assembly_v1.0.gff3 -o Acerv.merged stringtie_acerv_merged.gtf
48478 reference transcripts loaded.
  101293 query transfrags loaded
```

Now look at the summary file (Acerv.merged): 

<img width="799" alt="Screen Shot 2023-07-05 at 11 39 04 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/00aa7b4b-c412-4dd2-bcf3-47712bf62ff9">
