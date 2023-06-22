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

Ok now let's try the most recent alignment code but with the locally installed STAR.

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

for sample in ${data[@]} ;

do \
/scratch/projects/and_transcriptomics/programs/STAR-2.7.10b/bin/Linux_x86_64_static/STAR \
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

Ok this time when I submitted it, I think it couldn't read the fastq files because they were still in the zipped .gz format (which is why in Natalia's script she has a flag for --readFilesCommand gunzip -c). I'll try adding that then running it again.

Ok I think that's working.

Update: It did work!!!! I did get an exit code though because I think there was an issue with the multiqc part of the code. I removed the last two lines of the code to get rid of the multiqc code, and I'll try running that as a separate script now.

This is what I ran and it worked:
```{bash}
#!/bin/bash
#BSUB -J AcerCCC_multiqc
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o multiqc_STARalign.out
#BSUB -e multiqc_STARalign.err
#BSUB -n 8
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

cd /scratch/projects/and_transcriptomics/Allyson_CCC/aligned/

multiqc .
```

I can't tell if the results of the alignment are good or not.

![Screen Shot 2023-06-22 at 1 09 22 PM](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/43e0e9b0-e61d-41fe-a53f-fe0fdae85e5e)


![Screen Shot 2023-06-22 at 1 09 36 PM](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/bbc35f87-52b3-4104-b71a-83f8f1a2db48)


![Screen Shot 2023-06-22 at 1 09 48 PM](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/59a8a10e-1349-483f-95a1-43fb5df772f5)

If I had to guess, I would say the results are not good. 

Thankfully, Jill Ashey has come to save the day by sharing her github repo for using STAR + the Acer genome from Iliana Baums and Sheila Kitchen. https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/RNASeq_pipeline_FL.md 

In her code, one thing that could be useful and help is this note that she added: 
![Screen Shot 2023-06-22 at 1 13 00 PM](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/08eeb8d7-4609-4908-b414-425ce3ffc197)

The only problem is I'm not sure how exactly to add the identifier. But she definitely added it for Acer before running STAR.

Here is the code she used for the Acer genome index:
```{bash}
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Acerv --genomeFastaFiles /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0_171209.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Acerv/Acerv.GFFannotations.fixed_transcript.gff3
```

and for aligning reads to the genome
```{bash}
mkdir AlignReads_Acerv
cd AlignReads_Acerv
ln -s /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/*trim.fq .

nano AlignReads_acerv.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="Align_Acerv_out_error"
#SBATCH --output="Align_Acerv_out"

module load STAR/2.5.3a-foss-2016b

F=/data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Acerv

array1=($(ls $F/*trim.fq))
for i in ${array1[@]}
do
STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir ${i}_TMP --readFilesIn ${i} --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Acerv --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix ${i}.
done 

sbatch AlignReads_acerv.sh 
```

