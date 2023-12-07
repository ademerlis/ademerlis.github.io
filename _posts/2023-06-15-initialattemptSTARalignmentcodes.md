---
layout: post
title: initial attempt STAR alignment codes
date: '2023-06-15'
categories: coding
tags: [coding, CCC_ch4, temperaturevariability2023]
---

Here is the script Natalia used to set up the STAR genome index (https://github.com/China2302/SCTLD_RRC/blob/main/hpc/STAR_index.sh):
```{bash}
#!/bin/bash
#BSUB -J star_index
#BSUB -q general
#BSUB -P c_transcriptomics
#BSUB -n 8
#BSUB -o /scratch/projects/c_transcriptomics/Genomes/Ofav/star_index%J.out
#BSUB -e /scratch/projects/c_transcriptomics/Genomes/Ofav/star_index%J.err
#BSUB -u nxa945@miami.edu
#BSUB -N

todata="/scratch/projects/c_transcriptomics"
STAR \
--runThreadN 16 NumberOfThreads \
--runMode genomeGenerate \
--genomeDir ${todata}/Genomes/Ofav/Ofav_index \
--genomeFastaFiles ${todata}/Genomes/Ofav/GCF_002042975.1_ofav_dov_v1_genomic.fna \
--sjdbGTFfile /${todata}/Genomes/Ofav/genomic.gtf \
--sjdbOverhang ReadLength-100 \
--genomeSAindexNbases 13
```

First, I need to install STAR? or check that it exists on Pegasus.

Yes, it does exist on Pegaus. /share/Modules/bio/star/2.3.0e

So, to load it, follow the same code you used for fastqc (because this also exists on pegasus already)

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
--runThreadN 16 NumberOfThreads \
--runMode genomeGenerate \
--genomeDir ${and}/genomes/Acer/Acer_STAR_index \
--genomeFastaFiles ${and}/genomes/Acer/Acerv_assembly_v1.0_171209.fasta \
--sjdbGTFtagExonParentTranscript /${and}/genomes/Acer/Acerv_assembly_v1.0.gff3 \
--sjdbOverhang ReadLength-100 \
--genomeSAindexNbases 13
```

I'm going to try this and see if it works.

So far it's running and not immediately quitting.

Following the index generation, we then run a script to align reads to the indexed genome.

Here's Natalia's code (https://github.com/China2302/SCTLD_RRC/blob/main/hpc/star_align_trimmed.sh):
```{bash}
#!/bin/bash
#BSUB -J star_align
#BSUB -q bigmem
#BSUB -P c_transcriptomics
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -o star_align%J.out
#BSUB -e star_align%J.err
#BSUB -u nxa945@miami.edu
#BSUB -N

# Data used for this allignment was trimmed with "timming.sh" parameters. 
# Notice that the sequences still have a lot of polyA and adapter contamination
# a soft cliping option is added to STAR to deal with it.

todata="/scratch/projects/c_transcriptomics/"

cd "/nethome/nxa945/forStar"

data=($(ls *.gz))

for sample in ${data[@]} ;

do \
STAR \
--runMode alignReads \
--genomeDir ${todata}/Genomes/Ofav/Ofav_index \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfile  ${todata}/Genomes/Ofav/genomic.gtf \
--twopassMode Basic \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 16 \
--readFilesCommand gunzip -c \
--readFilesIn ${sample} \
--outFilterMultimapNmax 20 \
--quantMode TranscriptomeSAM GeneCounts \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--clip3pAdapterSeq AAAAAAAAAA \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outReadsUnmapped Fastx \
--outFileNamePrefix ${todata}/NOVA_SCTLD/data/alligned/${sample} ; \

done

multiqc ${todata}/NOVA_SCTLD/data/alligned/ \
--outdir ${todata}/NOVA_SCTLD/data/alligned/
```

For mine:
```{bash}
#!/bin/bash
#BSUB -J star_align_trimmed
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

and="/scratch/projects/and_transcriptomics/"

cd "/scratch/projects/and_transcriptomics/Allyson_CCC/trimmed/trimmed_and_removedpolyA_fastqfiles/forSTAR"

data=($(ls *.gz))

for sample in ${data[@]} ;

do \
module load star/2.3.0e \
STAR \
--runMode alignReads \
--genomeDir ${and}/genomes/Acer/Acer_STAR_index \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript Parent \
--twopassMode Basic \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 16 \
--readFilesCommand gunzip -c \
--readFilesIn ${sample} \
--outFilterMultimapNmax 20 \
--quantMode TranscriptomeSAM GeneCounts \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--clip3pAdapterSeq AAAAAAAAAA \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outReadsUnmapped Fastx \
--outFileNamePrefix ${and}/Allyson_CCC/aligned/${sample} ; \

done

multiqc ${and}/Allyson_CCC/aligned/ \
--outdir ${todata}/Allyson_CCC/aligned/
```

Let's see if this works.

Okay that didn't work. I think there is an issue whenever I try to load a module within the for-do loop. Trying this now:

```{bash}
#!/bin/bash
#BSUB -J star_align_trimmed
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

module load star/2.3.0e

for sample in ${data[@]} ;

do \
STAR \
--runMode alignReads \
--genomeDir ${and}/genomes/Acer/Acer_STAR_index \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript Parent \
--twopassMode Basic \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 16 \
--readFilesCommand gunzip -c \
--readFilesIn ${sample} \
--outFilterMultimapNmax 20 \
--quantMode TranscriptomeSAM GeneCounts \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--clip3pAdapterSeq AAAAAAAAAA \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outReadsUnmapped Fastx \
--outFileNamePrefix ${and}/Allyson_CCC/aligned/${sample} ; \

done

multiqc ${and}/Allyson_CCC/aligned/ \
--outdir ${todata}/Allyson_CCC/aligned/
```

This time i got this error: "FATAL INPUT ERROR: unrecoginzed parameter name "twopassMode""

I removed "--twopassMode Basic \" to see if that worked, and then I got another error based on a different flag. If I had to guess, I would say it is either an issue with the version of STAR that's on Pegasus is outdated (the author Dobin recommended to compile program from source as it has the most recent debugs that would fix this issue: https://groups.google.com/g/rna-star/c/-yWY0YB3dSM). Or, it could be that because I'm not providing an annotation file (.gtf file), none of these things work. What is frustrating is that the Acer genome files I downloaded from Galaxy don't have a gtf file (https://usegalaxy.org/published/history?id=1f8678b27ae56467), just a txt file for annotations. I looked into whether that could be directly converted and it isn't clear. 

But then looking at the older Acer genome that's available online (https://www.dropbox.com/s/wovxoi2dxcz9kv2/a.cervicornis_2014.zip?dl=0&file_subpath=%2Fa.cervicornis) it also doesn't have a .gtf file. 

I guess I'll try to download the latest version of STAR first? Or look for other Acer transcriptomic pipelines.
