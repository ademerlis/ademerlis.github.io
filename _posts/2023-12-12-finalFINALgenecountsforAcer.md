---
layout: post
title: Acer gene counts using samtools
date: '2023-12-07'
categories: Coding
tags:
  - Coding
published: true
---

## Finally getting the gene counts for Acer samples from Chapter 2

Following the last section of [Michael Studivan's code](https://github.com/mstudiva/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt). 

```{bash}
# NOTE: Must have a tab-delimited file giving correspondence between contigs in the transcriptome fasta file and genes

cat Acer/Locatelli_2023/Acervicornis_seq2iso.tab Symbiodinium/Mstudiva_Symbiodinium_annotated_transcriptome_code/Symbiodinium_seq2iso.tab > Host_concat_seq2iso.tab
```

Job: use sam tools and [samcount.pl](https://github.com/mstudiva/tag-based_RNAseq/blob/master/samcount.pl).

arg1: SAM file (by cluster, contig, or isotig)
arg2: a table in the form 'reference_seq<tab>gene_ID', giving the correspondence of 
reference sequences to genes. With 454-deived transcriptome, the gene_ID would be isogroup; 
with Trinity-derived transcriptiome,it would be component.
  
paths:
  1. /scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer/bowtie2align_LocatelliShoguchi/sam_files/*.sam
  2. /scratch/projects/and_transcriptomics/genomes/Host_concat_seq2iso.tab
  
I need to write a for-loop to each .sam file to do this: print "samcount.pl $file $tab aligner=bowtie2 >$file.counts" but written in bash not perl.
  
```{bash}
#! /usr/bin/env bash

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir=

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer/bowtie2align_LocatelliShoguchi/sam_files"

data=($(ls *.sam))

for samp in "${data[@]}" ; do \

#build script
echo "making sam_counts script for ${samp}..."
echo "
#! /usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${samp}_samcounts
#BSUB -e ${and}/Ch2_temperaturevariability2023/3_genecounts/logs/${samp}_samcounts.err
#BSUB -o ${and}/Ch2_temperaturevariability2023/3_genecounts/logs/${samp}_samcounts.out
#BSUB -W 12:00
#BSUB -n 8
#BSUB -q general

cd \"/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer/bowtie2align_LocatelliShoguchi/sam_files\"

module load samtools/1.3

perl samcount.pl ${samp} /scratch/projects/and_transcriptomics/genomes/Host_concat_seq2iso.tab aligner=bowtie2 >${samp}.counts

" > ${and}/Ch2_temperaturevariability2023/3_genecounts/${samp}_samcounts.job

bsub < ${and}/Ch2_temperaturevariability2023/3_genecounts/${samp}_samcounts.job

done
```

When I ran this though on Pegasus, the counts files came out with nothing in them. Moreover, all the .err files have errors like:

<img width="959" alt="Screen Shot 2023-12-15 at 11 33 09 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/af05892b-22cd-4265-8d3f-d37aa7cd5e33">

Why didn't any of this work? It looks like it maybe has to do with the isogroup assignment? 

This is what Michael ran on his samples before making the seq2iso.tab file: ([link](https://github.com/mstudiva/Acropora-cervicornis-annotated-transcriptome/blob/main/tagSeq_TranscriptomeAnnotation_README.txt))

```{bash}
# Renaming gene identifiers for ease
sed -i 's/comp/Acropora/g' Acervicornis.fasta
sed -i 's/isogroup/Acropora/g' Acervicornis.fasta
```

I didn't initially run this on my fasta file but I'll try it now and see if it helps.

```{bash}
# Renaming gene identifiers for ease
sed -i 's/comp/Acropora/g' Acropora_cervicornis.mrna-transcripts.fa > Acer2023_edited.fasta
sed -i 's/isogroup/Acropora/g' Acropora_cervicornis.mrna-transcripts.fa > Acer2023_edited.fasta
```

It doesn't help because there is nothing labeled "comp" or "isogroup in the fasta file I'm using from Locatelli, which is called Acropora_cervicornis.mrna-transcripts.fa. I looked through the other file types (.gff3, protein.fa, scaffolds.fa) and none of them have the comp or isogroup headers.

I ended up running this code to create a seq2iso.tab file.

```{bash}
grep ">" Acropora_cervicornis.mrna-transcripts.fa | perl -pe 's/>FUN(\d+)(\S+)\s.+/FUN$1$2\tFUN$1/'>Acervicornis_seq2iso.tab
```
This essentially just extracts the header line from each gene and puts it in a tab-delimited text file. 

<img width="384" alt="Screen Shot 2023-12-20 at 9 04 34 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/a07fbb1b-cc5b-4361-a87c-0d2949be5f32">

But then when I run samcount.pl using this seq2iso.tab file, I get no results.

**December 20 update:**

I had Michael take a look at the seq2iso.tab file, and the issue is that each line starts with ">" but it shouldn't. It's messing up the gene names.

So two things I need to do: First, update the Acropora one to have the genes be named "Acropora_000001" instead of "FUN_000001", and then update the Symbiodinium fasta file to say Symbiodinium_000001 for those genes. Then create the seq2iso.tab files for both species, concatenate them, then re-run the samcount.pl job.

```{bash}
# Renaming gene identifiers for ease

#make new file
cp Acropora_cervicornis.mrna-transcripts.fa Acer_edited.fasta
cp syma_transcriptome_37.fasta Symbiodinium_edited.fasta

sed -i 's/FUN/Acropora/g' Acer_edited.fasta

sed -i 's/comp/Symbiodinium/g' Symbiodinium_edited.fasta
sed -i 's/isogroup/Symbiodinium/g' Symbiodinium_edited.fasta


# making seq2iso.tab files
grep ">" Acer_edited.fasta | perl -pe 's/>Acropora_(\d+)(\S+)\s.+/Acropora_$1$2\t Acropora_$1/'>Acervicornis_seq2iso.tab
grep ">" Symbiodinium_edited.fasta | perl -pe 's/>Symbiodinium(\d+)(\S+)\s.+/Symbiodinium$1$2\t Symbiodinium$1/'>Symbiodinium_seq2iso.tab

# create combo file

cat Acer/Locatelli_2023/Acer_Genome/Acervicornis_seq2iso.tab Symbiodinium/Symbiodinium_seq2iso.tab > Host_concat_seq2iso.tab
```

~~Now i should be able to re-run this script with no changes:~~

Incorrect, the sam files are all named with the original gene names so nothing is going to match up.

I can do one of two things: either rename all the same files to have Acropora instead of FUN, or remake the seq2iso.tab files with the original names and then rename everything after the fact when its in the gene counts matrix. 

The sam files are huge so I'm going to try to rename them after the fact. 

```{bash}
# making seq2iso.tab files
grep ">" Acropora_cervicornis.mrna-transcripts.fa | perl -pe 's/>FUN_(\d+)(\S+)\s.+/FUN_$1$2\t FUN_$1/'>Acervicornis_seq2iso.tab
grep ">" syma_transcriptome_37.fasta | perl -pe 's/>comp(\d+)(\S+)\s.+/comp$1$2\t comp$1/'>Symbiodinium_seq2iso.tab

# create combo file

cat Acer/Locatelli_2023/Acer_Genome/Acervicornis_seq2iso.tab Symbiodinium/Symbiodinium_seq2iso.tab > Host_concat_seq2iso.tab
```
