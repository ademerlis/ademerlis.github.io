---
layout: post
title: HISAT2 alignment
date: '2023-07-23'
categories: coding
tags: [coding, Ch2_tempvariability]
---

So, following the next step of the pipeline, Ariana, Kevin, Sam, and Zoe all use HISAT2 for alignment. Here is an explanation of all the flags from [Sam's GitHub](https://github.com/SamGurr/SamGurr.github.io/blob/master/_posts/2021-01-07-Geoduck-TagSeq-Pipeline.md#hisat2-alignment). 

Samtools is a module on Pegasus, but HISAT2 is not. Both are called in the job scripts so I need to install HISAT2 to my scratch space.

I installed it from [source](http://daehwankimlab.github.io/hisat2/download/#version-hisat2-221) (locally downloaded then scp onto Pegasus and added PATH to ~/.bash_profile). Also need to run "make" in the terminal. (this creates important executables like hisat-align-1)

HISAT2 works just like STAR. First, you need to index the genome, then you align reads to the reference genome.

Sam does both steps in one code:
```{bash}
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=200GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=samuel_gurr@uri.edu
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"
#SBATCH -D /data/putnamlab/sgurr/Geoduck_TagSeq/output/HISAT2

#load packages
module load HISAT2/2.1.0-foss-2018b #Alignment to reference genome: HISAT2
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

# symbolically link 'clean' reads to hisat2 dir
ln -s ../fastp_multiQC/clean*.fastq.gz ./ 

# index the reference genome for Panopea generosa output index to working directory
hisat2-build -f ../../../refs/Panopea-generosa-v1.0.fa ./Pgenerosa_ref # called the reference genome (scaffolds)
echo "Referece genome indexed. Starting alingment" $(date)

# This script exports alignments as bam files
# sorts the bam file because Stringtie takes a sorted file for input (--dta)
# removes the sam file because it is no longer needed
array=($(ls *.fastq.gz)) # call the symbolically linked sequences - make an array to align
for i in ${array[@]}; do
        sample_name=`echo $i| awk -F [.] '{print $2}'`
	hisat2 -p 8 --dta -x Pgenerosa_ref -U ${i} -S ${sample_name}.sam
        samtools sort -@ 8 -o ${sample_name}.bam ${sample_name}.sam
    		echo "${i} bam-ified!"
        rm ${sample_name}.sam
done
```

I will adapt the above code to run on Pegasus:
```{bash}
#!/bin/bash
#BSUB -J HISAT2
#BSUB -q bigmem
#BSUB -n 16
#BSUB -P and_transcriptomics
#BSUB -o HISAT2%J.out
#BSUB -e HISAT2%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

module load samtools/1.3
module load python/3.8.7

/scratch/projects/and_transcriptomics/programs/hisat2-2.2.1/hisat2-build -f ${and}/genomes/Acer/Acerv_assembly_v1.0_171209.fasta ${and}/genomes/Acer/Acer_STAR_index_gffannotations.fixed_take3
echo "Reference genome indexed. Starting alignment" $(date)

array=($(ls /scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/AS_pipeline/trimmed/*.fastq.gz))
for i in ${array[@]}; do
        sample_name=`echo $i| awk -F [.] '{print $2}'`
	/scratch/projects/and_transcriptomics/programs/hisat2-2.2.1/hisat2 -p 8 --dta -x ${and}/genomes/Acer/Acer_STAR_index_gffannotations.fixed_take3 -U ${i} -S ${sample_name}.sam
        samtools sort -@ 8 -o ${sample_name}.bam ${sample_name}.sam
    		echo "${i} bam-ified!"
        rm ${sample_name}.sam
done
```

It's weird, I keep getting an error that the trimmed files (i.e. "Acer-005_S23_L001_R1_001.fastq.gz") are not in gzipped format. But they have the end of ".gz". It makes me think I messed up the files somehow when I did the trimming. I'm going to re-run that step and make sure the naming convention doesn't add modifiers to the end of the file names.


