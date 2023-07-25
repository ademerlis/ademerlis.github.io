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

I re-ran the trimming step and HISAT2 was able to use those files so that's good.

But then I also got this error: "(ERR): hisat2-align exited with value 1
/scratch/projects/and_transcriptomics/programs/hisat2-2.2.1/hisat2-align-s: invalid option -- '@'

I think what fixed it is I changed the script so there weren't any erroneous "\" symbols -- I think when I use that, it means to continue on with the second line of the script as if it was with the first line. But we want it to do the functions sequentially, not at the same time. So the -- '@' flag should correspond to samtools, not to hisat2-align. 

So I ran it again and I think I was getting it to work (see below screenshot for .err file):

<img width="311" alt="Screen Shot 2023-07-25 at 9 11 49 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/55af7ccc-c3d2-4221-96dc-67ab6fd44c2e">

LOL it tried to align the Pcli samples too. It got 0.01% alignment lol. 

To view mapping percentages (percent alignment): ([code from Zoe](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/ZD_Heron-Pdam-gene-expression.md))

```{bash}
module load samtools/1.3

for i in *.bam; do
    echo "${i}" >> mapped_reads_counts_Acer
    samtools flagstat ${i} | grep "mapped (" >> mapped_reads_counts_Acer
done
```

The output looks like this, but I don't know what the other numbers mean (besides the percent alignment).

<img width="218" alt="Screen Shot 2023-07-25 at 9 35 29 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/06ec87c2-e808-41b1-adfa-0fbecd0c9929">

Ok when I run "samtools flagstat" on one file, I get all this info, but none of it looks super useful except for what was already extracted in the for loop above:

<img width="551" alt="Screen Shot 2023-07-25 at 9 37 35 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/f69d4ae4-d7a6-47c9-ab6e-bed58c4895c4">

I think the number that I was confused about was the X + 0 mapped, and that is separating QC-passed reads + QC-failed reads. 

cool so I think the alignment went well! Time for stringtie.

But before that, I just wanted to compare alignment rates of STAR versus HISAT2 and see how different it was. 

What's annoying is that from hisat2 you can't multiqc it like you can from STAR (strike against hisat2 in my mind). 

So here's a rough screenshot of the HISAT2 results in the terminal on the left, and the STAR alignment results on the right (for a subset of samples):

<img width="926" alt="Screen Shot 2023-07-25 at 9 42 01 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/15fb67eb-1949-4dec-a9b4-628f407b2618">

Overall alignment rates look pretty similar so that's good. 
