---
layout: post
title: stringtie code for Acer CCC
date: '2023-06-27'
categories: coding
tags: [coding, ch4_CCC]
---

1) Need to install stringtie and gffcompare to local programs folder on pegasus.

```{bash}
wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.1.tar.gz
tar xvfz stringtie-2.2.1.tar.gz
cd stringtie-2.2.1
make release

wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.12.6.Linux_x86_64.tar.gz
tar xvfz gffcompare-0.12.6.Linux_x86_64.tar.gz
```
For stringtie i downloaded the source version then had to "make" it, but i could've also downloaded the pre-compiled version which is what i did for gffcompare (the Linux_x86_64 version). So i don't need to "make" or "source' it.

2) Ok now that I have both downloaded I can start making my script. NOTE: all these scripts are adapted from the wonderful [Jill Ashey](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/RNASeq_pipeline_FL.md)

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

cd "/scratch/projects/and_transcriptomics/Allyson_CCC/aligned"

data=($(ls *Aligned.sortedByCoord.out.bam))

for i in ${data[@]} ;

do \
${and}/programs/stringtie-2.2.1/stringtie -G /${and}/genomes/Acer/Acerv_assembly_v1.0.gff3 -e -o ${i}.gtf ${i} ; \ 
done
```

Ok what do all the flags mean:
-G is a reference annotation file to be used as a guide for the assembly process. 
-e "When the -e option is used, the reference annotation file -G is a required input and StringTie will not attempt to assemble the input read alignments but instead it will only estimate the expression levels of the "reference" transcripts provided in the -G file. With this option, no "novel" transcript assemblies (isoforms) will be produced, and read alignments not overlapping any of the given reference transcripts will be ignored, which may provide a considerable speed boost when the given set of reference transcripts is limited to a set of target genes for example."
-o output file name

I think this is working (although it says "line 21: command not found") but there are still .gtf files being generated so that's good I guess.

3) The next step is then to **Merge stringTie gtf results**:

```{bash}
mkdir stringtie_gtf_files
mv *.gtf stringtie_gtf_files/
ls *gtf > acerv_mergelist.txt
cat acerv_mergelist.txt
```

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

cd "/scratch/projects/and_transcriptomics/Allyson_CCC/aligned/stringtie_gtf_files"

${and}/programs/stringtie-2.2.1/stringtie --merge -p 8 -G ${and}/genomes/Acer/Acerv_assembly_v1.0.gff3 -o stringtie_acerv_merged.gtf acerv_mergelist.txt
```
That worked first try, nice

4) next is to **assess assembly quality** using gffcompare (i'm going to try not writing this as a job and just run it in the command line because these things seem to be running fast):

```{bash}
cd /scratch/projects/and_transcriptomics/Allyson_CCC/aligned/stringtie_gtf_files

/scratch/projects/and_transcriptomics/programs/gffcompare-0.12.6.Linux_x86_64/gffcompare -r /scratch/projects/and_transcriptomics/genomes/Acer/Acerv_assembly_v1.0.gff3 -o Acerv.merged stringtie_acerv_merged.gtf

 48478 reference transcripts loaded.
  48478 query transfrags loaded.
```

When I look at the summary of gffcompare, I see this (which strikes me as odd):

![Screen Shot 2023-06-27 at 3 50 40 PM](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/daf42997-88e7-4876-9e5b-d56f678a9a1e)

I don't know what the numbers should be, but i feel like having 100% for everything means it's not reading the data correctly. 

Things that could be why it didn't work:
- maybe i need to include the full path to each of the .bam.gtf files in the merged_list.txt there is a note in [Danielle's code - line 900](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-04-14-Molecular-Underpinnings-RNAseq-Workflow.md)
- i should also include the fasta file [pg 14 of this doc](https://bioinformatics-core-shared-training.github.io/RNAseq_September_2018/slides/DAY_1.pdf)
- I found this [github issue thread](https://github.com/gpertea/gffcompare/issues/21) that might be the same problem I'm having. It goes back to the "transcript_id" thing that Jill had to correct in the Acer genome assembly .gff3 file before runnning STAR. Maybe because I didn't do that, I'm getting 100% alignment because everything is under the same column or something.
- 
I found [Jill's code to fix the .gff3 file](https://github.com/JillAshey/SedimentStress/blob/master/RAnalysis/gff_wrangling/acerv_Fix_GFF.R), so I'm going to adapt that and run it in R. I think I can run an R environment on pegasus... The only thing is that I would need to install all the packages probably because idk what libarires Pegasus has installed. So I'll just download the .gff3 file locally, edit in R, then re-upload.  

```{r}
#from Jill
# Title: A. cervicornis GFF adjustments
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 09/01/20

# Need to do some acerv gff adjustments so it can run properly in STAR. Here, I'll be adding transcript_id= to 'gene' column because STAR needs that label to run

#Load libraries
library(tidyverse)

#Load  gene gff
Acerv.gff <- read.csv(file="~/Desktop/GFFs/Acerv_assembly_v1.0.gff3", header=FALSE, sep="\t", skip=1) 

#rename columns
colnames(Acerv.gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
head(Acerv.gff)

# Creating transcript id
Acerv.gff$transcript_id <- sub(";.*", "", Acerv.gff$gene)
Acerv.gff$transcript_id <- gsub("ID=", "", Acerv.gff$transcript_id) #remove ID= 

# Checking what kinds of ids are in gff
unique(Acerv.gff$id)
# [1] "gene"        "mRNA"        "exon"        "CDS"         "start_codon" "stop_codon"  "tRNA"       

#If id == mRNA, exon, start_codon, stop_codon, CDS, tRNA, add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with original gene
Acerv.gff <- Acerv.gff %>% 
  mutate(gene = ifelse(id != "gene", paste0(gene, ";transcript_id=", Acerv.gff$transcript_id),  paste0(gene)))
head(Acerv.gff)

# Remove last col
Acerv.gff <- Acerv.gff[,-10]
head(Acerv.gff)  

#save file
write.table(Acerv.gff, file="~/Desktop/GFFs/Acerv.GFFannotations.fixed_transcript.gff3", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)
```



5) **re-estimate assembly**
   (jill submitted this as a job so it must take some more computing power)
   
```{bash}
#!/bin/bash
#BSUB -J stringtie_reestimateassembly
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o stringtie_reestimateassembly%J.out
#BSUB -e stringtie_reestimateassembly%J.err
#BSUB -u and128@miami.edu
#BSUB -N

module load python/3.8.7

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Allyson_CCC/aligned"

data=($(ls *Aligned.sortedByCoord.out.bam))

for i in ${data[@]} ;

do \
${and}/programs/stringtie-2.2.1/stringtie -e -G /${and}/genomes/Acer/Acerv_assembly_v1.0.gff3 -o ${i}.merge.gtf ${i} ; \ 
done
```


