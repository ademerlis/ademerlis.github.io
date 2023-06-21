---
layout: post
title: STAR alignment Acer genome
date: '2023-06-21'
categories: coding
tags: [coding, CCC_ch4, temperaturevariability2023]
---

This was the most recent code I tried that didn't work:

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
```

All of these flags were copied from Natalia's code. But I'm not sure if all of these are needed?

The ones used for the Apal genome from Young et al. 2020 were:

```{bash}
echo '/nethome/bdy8/Ben_Xaymara_GE_project/programs/STAR \
--runThreadN 8 \
--genomeDir /scratch/projects/transcriptomics/ben_young/apalm_v2/star_index/ \
--readFilesIn /scratch/projects/transcriptomics/ben_young/apalm_v2/trimmed_files/'"$PALPAL"'_trimmed.fastq \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMstrandField intronMotif \
--twopassMode Basic \
--twopass1readsN -1 \
--outFileNamePrefix /scratch/projects/transcriptomics/ben_young/apalm_v2/aligned_apal/'"$PALPAL"'/'"$PALPAL"'_' >> /scratch/projects/transcriptomics/ben_young/apalm_v2/scripts/star_scripts/"$PALPAL"_alignment.sh #for alignment to symb
```

Maybe I should just try using the bare minimum flags and see if it works?

What's weird is "--runMode" isn't used in the second example.

I tried doing a frankenstein version of the two scripts and I still got an error (EXITING: FATAL INPUT ERROR: unrecoginzed parameter name "outSAMtype" in input "Command-Line-Initial"
SOLUTION: use correct parameter name (check the manual)).
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
--genomeDir ${and}/genomes/Acer/Acer_STAR_index/ \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 16 \
--readFilesIn ${sample} \
--quantMode TranscriptomeSAM GeneCounts \
--clip3pAdapterSeq AAAAAAAAAA \
--outReadsUnmapped Fastx \
--outFileNamePrefix ${and}/Allyson_CCC/aligned/${sample} ; \

done

multiqc ${and}/Allyson_CCC/aligned/ \
--outdir ${todata}/Allyson_CCC/aligned/
```

I'm just going to download the STAR program locally on my Pegasus scratch space and use the latest version. Maybe the one on Pegasus isn't the best version.

I should probably rerun the index script as well once I download the newest version of the STAR program so everything is compatible. 

## Installing STAR to scratch space ##

```{bash}
# Get latest STAR source from releases
#ran this all in login node 
wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz
tar -xzf 2.7.10b.tar.gz
cd STAR-2.7.10b
# Compile
cd source
make STAR
```

I think that worked. STAR is located at this absolute path: /scratch/projects/and_transcriptomics/programs/STAR-2.7.10b

I'll try to do the index again with this version now. 

Ok so it didn't work at first because the path to run STAR wasn't correct. I need to use this path when I want to load and run STAR locally within a script: /scratch/projects/and_transcriptomics/programs/STAR-2.7.10b/bin/Linux_x86_64_static/STAR

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

I got this warning: "--genomeSAindexNbases 14 is too large for the genome size=318373619, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 13"

So add that flag to the script and run it again.

Ok that seems to be working fine now. I'll download that script to my github files so I have the most up-to-date version.
