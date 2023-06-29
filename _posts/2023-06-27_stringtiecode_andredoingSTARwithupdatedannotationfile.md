---
layout: post
title: stringtie code for Acer CCC and rerunning STAR with updated annotation file
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
- G is a reference annotation file to be used as a guide for the assembly process. 
- e "When the -e option is used, the reference annotation file -G is a required input and StringTie will not attempt to assemble the input read alignments but instead it will only estimate the expression levels of the "reference" transcripts provided in the -G file. With this option, no "novel" transcript assemblies (isoforms) will be produced, and read alignments not overlapping any of the given reference transcripts will be ignored, which may provide a considerable speed boost when the given set of reference transcripts is limited to a set of target genes for example."
- o output file name

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

Because she did this step before the STAR alignment, I should re-run all the STAR things so that those files aren't messed up. 

**new STAR index**
```{bash}
#!/bin/bash
#BSUB -J Acer_star_index_fixedannotations
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -o /scratch/projects/and_transcriptomics/genomes/Acer/star_index%J.out
#BSUB -e /scratch/projects/and_transcriptomics/genomes/Acer/star_index%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

/scratch/projects/and_transcriptomics/programs/STAR-2.7.10b/bin/Linux_x86_64_static/STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeDir ${and}/genomes/Acer/Acer_STAR_index \
--genomeFastaFiles ${and}/genomes/Acer/Acerv_assembly_v1.0_171209.fasta \
--sjdbGTFfile ${and}/genomes/Acer/Acerv.GFFannotations.fixed_transcript.gff3 \
--sjdbOverhang 100 \
--sjdbGTFtagExonParentTranscript Parent \
--genomeSAindexNbases 13
```

**new STAR alignment**
```{bash}
#!/bin/bash
#BSUB -J star_align_trimmed_fixedannotations
#BSUB -q bigmem
#BSUB -P and_transcriptomics
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -o star_align%J.out
#BSUB -e star_align%J.err
#BSUB -u and128@miami.edu
#BSUB -N

# A soft clipping option is added to STAR to deal with any leftover polyA and adapter contamination 

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Allyson_CCC/trimmed/trimmed_and_removedpolyA_fastqfiles/forSTAR"

data=($(ls *.gz))

for sample in ${data[@]} ;

do \
/scratch/projects/and_transcriptomics/programs/STAR-2.7.10b/bin/Linux_x86_64_static/STAR \
--genomeDir ${and}/genomes/Acer/Acer_STAR_index_annotationsupdated \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 16 \
--readFilesCommand gunzip -c \
--readFilesIn ${sample} \
--quantMode TranscriptomeSAM GeneCounts \
--clip3pAdapterSeq AAAAAAAAAA \
--outReadsUnmapped Fastx \
--outFileNamePrefix ${and}/Allyson_CCC/aligned_updatedannotations/${sample} ; \

done
```

Then run multiqc again on the STAR alignment. 

**run stringtie again**
```{bash}
#!/bin/bash
#BSUB -J stringtie_updatedannotations
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o stringtie%J.out
#BSUB -e stringtie%J.err
#BSUB -u and128@miami.edu
#BSUB -N

module load python/3.8.7

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Allyson_CCC/aligned_updatedannotations"

data=($(ls *Aligned.sortedByCoord.out.bam))

for i in ${data[@]} ;

do \
${and}/programs/stringtie-2.2.1/stringtie -G /${and}/genomes/Acer/Acerv.GFFannotations.fixed_transcript.gff3 -e -o ${i}.gtf ${i} ; \ 
done
```

I still get the same error as  before (line 21: command not found) but the .gtf files are still generated so i'll move forward with them. 

**Merge stringTie gtf results**:

```{bash}
mkdir stringtie_gtf_files
mv *.gtf stringtie_gtf_files/
cd stringtie_gtf_files/
ls *gtf > acerv_mergelist.txt
cat acerv_mergelist.txt
```

```{bash}
#!/bin/bash
#BSUB -J stringtie_mergegtf_updatedannotations
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o stringtie_mergegtf_updatedannotations%J.out
#BSUB -e stringtie_mergegtf_updatedannotations%J.err
#BSUB -u and128@miami.edu
#BSUB -N

module load python/3.8.7

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Allyson_CCC/aligned_updatedannotations/stringtie_gtf_files"

${and}/programs/stringtie-2.2.1/stringtie --merge -p 8 -G ${and}/genomes/Acer/Acerv.GFFannotations.fixed_transcript.gff3 -o stringtie_acerv_merged.gtf acerv_mergelist.txt
```

**gffcompare to check assembly quality**

Can just run this in the terminal, it doesn't take long.
```{bash}
cd /scratch/projects/and_transcriptomics/Allyson_CCC/aligned_updatedannotations/stringtie_gtf_files

/scratch/projects/and_transcriptomics/programs/gffcompare-0.12.6.Linux_x86_64/gffcompare -r /scratch/projects/and_transcriptomics/genomes/Acer/Acerv.GFFannotations.fixed_transcript.gff3 -o Acerv.merged stringtie_acerv_merged.gtf

 48478 reference transcripts loaded.
  48478 query transfrags loaded.
```

Shoot I still get the same thing with 100% accuracy and precision...

<img width="522" alt="Screen Shot 2023-06-29 at 11 33 43 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/46f2eb54-07d0-4306-acd6-83d795f8dbd5">

I don't know what could be wrong with this. I asked ChatGPT for help and it gave me this: 

"If your `gffcompare` results show 100% sensitivity and precision, it typically suggests that there might be an issue with the input data or the `gffcompare` command itself. Here are a few potential reasons for this outcome:

1. **Incorrect input files:** Ensure that you have provided the correct input files to `gffcompare`. The reference annotation file (GTF/GFF) should be specified as the first input, followed by the assembled transcript file (GTF/GFF). Check if you have accidentally switched the input order.

2. **Mismatched or incompatible annotation versions:** `gffcompare` requires the reference annotation file and the assembled transcript file to be based on the same genome version and annotation guidelines. If there is a version mismatch or inconsistency, it can lead to inaccurate results. Double-check that both files correspond to the same reference genome and annotation version.

3. **Errors in the assembly process:** If the assembled transcript file generated by StringTie or any other assembler contains significant errors or artifacts, `gffcompare` might not be able to accurately compare it with the reference annotation. Review the assembly steps and parameters to ensure proper transcript reconstruction.

4. **Inadequate reference annotation:** If the reference annotation file used is incomplete or lacks comprehensive coverage of the transcriptome, it can lead to artificially high sensitivity and precision values. Assess the quality and completeness of the reference annotation to verify its adequacy for comparison.

5. **Other factors:** There could be other factors influencing the results, such as specific characteristics of the dataset, the presence of highly repetitive regions, or unusual transcript isoforms. Further investigation and additional analysis may be necessary to identify the cause of the unexpected outcome.

To troubleshoot the issue, you can:

- Double-check the command and input file formats to ensure they are correct.
- Verify that the input files are compatible and correspond to the same annotation version.
- Review the transcript assembly process and parameters for potential errors or inconsistencies.
- Inspect the reference annotation for completeness and suitability for comparison.
- Consider running `gffcompare` with different settings or parameters to see if the results change.

If the problem persists or you need further assistance, it would be helpful to provide more specific details about the command you used and the data you are working with."

I checked the order of input files and I'm pretty sure they're correct. I'll try to reorder them and see if that helps.

```{bash}
cd /scratch/projects/and_transcriptomics/Allyson_CCC/aligned_updatedannotations/stringtie_gtf_files

/scratch/projects/and_transcriptomics/programs/gffcompare-0.12.6.Linux_x86_64/gffcompare -r /scratch/projects/and_transcriptomics/genomes/Acer/Acerv.GFFannotations.fixed_transcript.gff3 stringtie_acerv_merged.gtf -o Acerv.merged

 48478 reference transcripts loaded.
  48478 query transfrags loaded.
```
Nope that did nothing. 

The other thing I'm wondering is whether the merge function worked? Because the size of the merged gtf file is 44 MB while the size of each gtf file is ~42 MB. Shouldn't the merged one be way larger? But if the merged file only had one sample gtf in it, that still wouldn't make the senstivity and precision 100%... 

The multiqc for the updated STAR alignment came back as well and it's identical to the unedited gff one (in terms of alignment rates for the samples). So idk if that really did anything.

What if the gffcompare thing isn't working because of the error message I got? line 21: command not found? Maybe that's why it isn't working and it's generating the gtf files wrong or something. 

So I asked ChatGPT and apparently after Sept 2021 the -e option got removed from stringtie. So what does that mean the command is doing then? 

I just looked on the github page for stringtie and it says that the -e option is still there, but is required with -G. And the purpose of it is to "only estimate the abundance of given reference transcripts." Which idk what that means, like would that mean that the sample gtf files would only be matched to reference transcripts from the Baums and Kitchen genome assembly? 

Let me try two things: first try running it without the -e option and see if that does anything, then try running it with the fasta file as well to see if it needs like extra annotations or something. 

1. without e option
- even when I removed the e option it says "/projects/lsf_spool/1688064612.27955592.shell: line 21:  : command not found". So something isn't working at all. Maybe. I should try just runnnig stringtie on one file.
```{bash}
/scratch/projects/and_transcriptomics/programs/stringtie-2.2.1/stringtie -G /scratch/projects/and_transcriptomics/genomes/Acer/Acerv_assembly_v1.0.gff3 -o 1087_trimmed_trimmed.fastq.gzAligned.sortedByCoord.out.bam.gtf 1087_trimmed_trimmed.fastq.gzAligned.sortedByCoord.out.bam
```
It ran and didn't have any errors come up. Btu the resulting .gtf file is only 15KB in size which seems weird.

Let's try running it with the fasta file too and see if the file size changes. 

```{bash}
/scratch/projects/and_transcriptomics/programs/stringtie-2.2.1/stringtie -G /scratch/projects/and_transcriptomics/genomes/Acer/Acerv_assembly_v1.0.gff3 --ref /scratch/projects/and_transcriptomics/genomes/Acer/Acerv_assembly_v1.0_171209.fasta -o 1087_trimmed_trimmed.fastq.gzAligned.sortedByCoord.out.bam.gtf 1087_trimmed_trimmed.fastq.gzAligned.sortedByCoord.out.bam
```

It doesn't look like this made any difference. I guess I'll try to add this fasta reference thing into the script and see if the error message goes away. 

The error message didn't go away. Then I removed the .gff3 assembly file too and the ERROR STILL POPPED UP. so it must have something to do with stringtie program???

I just updated the ~/.bash_profile to have an export PATH specifically for stringtie. I can now run stringtie in the command line. But when I submit the job, i still get the error of line 21 command not found!!! what the heck am i doing wrong. 

I think I. amgoing to give up and just try Natalia's method of directly using the STAR output files. I'm still not sure what Stringtie even does to be honest.

 

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


