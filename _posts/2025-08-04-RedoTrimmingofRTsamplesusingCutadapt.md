---
layout: post
title: Redoing trimming of reciprocal transplant seqs
date: '2025-08-04'
categories: Analysis, Processing
tags: [Ch3, reciprocaltransplant, RNAseq, trimming, cutadapt, Pegasus]
---

I realized that I didn't follow the Lexogen guidelines for trimming sequences from samples that were prepped for cDNA libraries using QuantSeq. 

### Notes (June 19, 2025)
___ 

- Michael added me to a Box folder from UM Genomics center back in December 2024 - it is called “ReciprocalTransplant_11122024”
- each sample has two runs (R1=left and R2=right) which indicate that it’s paired-end reads.
- I guess in March 2025 I downloaded the sequences already onto the “2TB” external hard drive and uploaded them to Pegasus. (https://github.com/ademerlis/ademerlis.github.io/blob/main/_posts/2025-03-07-UploadingReciprocalTransplantSequencestoPegasus.md)
- So my 3 different dissertation chapters were sequenced differently, so they need to be analyzed differently.
    - Ch 2 - Acer and Pcli stress-hardening; sequencing done at UT Austin via 3’ Tag-Seq, so I used Michael’s cutadapt trimming specific code to make sure I removed the headers correctly; single-end reads
    - Ch 3 - Acer; sequencing done at UM Genomics Center through Natalia and Lys’s plate, QuantSeq (3’ RNA-Seq?) — so I used Natalia’s TrimGalore code to ensure I trimmed it how she did it; single-end reads
    - Ch 4 - Pstr; sequencing done at UM Genomics Center; paired-end reads
- How do I know what to do for these Pstr reads?
- Notes on sequencing from emails from Michael and Mia:
    - submit pooled RNA-Seq libraries for sequencing on 1 NovaSeq X PE150 lane per NOAA contract 1333MJ24P0082. Based on the pooling calculations, putting in 50 fmol per 96 samples will yield ~550 uL. 50 fmol per sample would result in a final pooled molarity of 8.72 nM in 550 uL.
    - the libraries were pooled and sequenced on an Illumina platform (Illumina NovaSeq X, PE150 mode)
- Notes from Lexogen (because cDNA libraries were prepared using QuantSeq FWD kit)
    - https://faqs.lexogen.com/faq/can-i-use-paired-end-sequencing-for-quantseq-fwd-l
        - “Paired-end (PE) sequencing is typically not recommended for QuantSeq FWD (Cat. No. 191 - 196), as the quality of Read 2 is very low due to the poly(T) stretch at the beginning of Read 2.
        - Nevertheless, QuantSeq FWD libraries can be sequenced in PE mode (i.e., 150bp) but in this case, we recommend discarding Read 2 data and proceeding with Read 1 data only for downstream data analysis (i.e., use only Read 1 for trimming, alignment, read counting, and downstream analyses).”
    - https://faqs.lexogen.com/faq/what-sequences-should-be-trimmed
        - they have specific cutadapt flags that they recommend

### Notes from ChatGPT (August 4, 2025)
___

So my sequences were prepared using the cDNA libary kit from Lexogen (**QuantSeq FWD kit**), and were sequenced on an Illumina platform (**NovaSeq X, PE150 mode**). 

ChatGPT said that there are trimming steps required for both QuantSeq and NovaSeq, which were not specifically captured in the TrimGalore script that I previously used (see below):

```{bash}
#! /usr/bin/env bash

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="/scratch/projects/and_transcriptomics/reciprocaltransplant"

cd "/scratch/projects/and_transcriptomics/reciprocaltransplant/raw_seq_files"

data=($(ls *.fastq.gz))

for samp in "${data[@]}" ; do \

#build script
echo "making TrimGalore script for ${samp}..."
echo "
#! /usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${samp}_trim
#BSUB -e ${projdir}/logs/${samp}_trim.err
#BSUB -o ${projdir}/logs/${samp}_trim.out
#BSUB -q general

cd \"/scratch/projects/and_transcriptomics/reciprocaltransplant/raw_seq_files\"

${and}/programs/TrimGalore-0.6.10/trim_galore ${samp} --illumina --cores 4 --three_prime_clip_R1 12 --three_prime_clip_R2 12 --nextseq 30 --length 20


" > ${projdir}/${samp}_trim.job

#submit script

bsub < ${projdir}/${samp}_trim.job

done
```

The flags `--illumina --cores 4 --three_prime_clip_R1 12 --three_prime_clip_R2 12 --nextseq 30 --length 20`are simple for trimming, but have a couple of limitations.

- The --illumina flag tells TrimGalore to look for Illumina adapters, but isn't specific.
- the --three_prime_clip hard-trims 12 bases from the 3' ends, but doesn't do a soft-trim.
- --nextseq 30 trims G tails due to NextSeq chemistry, which creates a poly-G artifact from 2-color sequencing, using a Q30 cutoff. This is important for Illumina NextSeq/NovaSeq
- --length just discards reads shorter than 20 bp.

Lexogen's trimming strategy does multiple passes, and is more specific with adapters.

Also, as per Lexogen's suggestion, read 2 (R2) is pretty much useless for downstream analysis, so I am not going to work on trimming those anymore. I'm going to move them all into a subfolder in Pegasus and stop processing them. 
- But this begs the question then, do people use them in de novo transcriptome sequencing?? When they're using paired-end reads, what kind of reads/sequencing were they doing to be able to do that?

I may need to complete my literature review before I keep going through this.

The lower-quality reads of R2 is in agreement with what is seen in the multiQC report even following trimming (see screenshot below and MultiQC report [here](https://github.com/ademerlis/reciprocaltransplant/blob/main/gene_expression/multiqc_reports/multiqc_report_trimmedreads.html)), at least for the mean quality scores. The 97 green sequences are all the R1s.

<img width="1453" height="668" alt="Screenshot 2025-08-04 at 4 10 57 PM" src="https://github.com/user-attachments/assets/a9bc84fe-225d-4788-92eb-c67400341fcd" />

