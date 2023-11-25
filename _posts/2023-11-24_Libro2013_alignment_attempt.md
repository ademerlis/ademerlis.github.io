---
layout: post
title: Libro et al. 2013 Acer genome alignment attempt
date: '2023-11-24'
categories: Coding
tags: [Coding]
---

I'm going to try to run the bowtie2 scripts from [Dr. Michael Studivan](https://github.com/mstudiva/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt) and see how the alignment rates come out. I'm going to concatenate the Acer genome with the Shoguchi et la. 2018 Symbiodinium genome. 

First, run bowtie2-build:

```{bash}
#!/bin/bash
#BSUB -J bowtie2-build_libro_shoguchi
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -W 120:00
#BSUB -o bowtie2-build_libro_shoguchi.out
#BSUB -e bowtie2-build_libro_shoguchi.err
#BSUB -u and128@miami.edu
#BSUB -N

workdir="/scratch/projects/and_transcriptomics/genomes"

bowtie2-build ${workdir}/Acer_2023/Libro_2013/acer.fasta, \
${workdir}/Symbiodinium/syma_transcriptome_37.fasta \
Host_concat
```

Next, run alignment.

Updated bowtie2 code based on ChatGPT's recommendations (some of the flags Michael was using aren't recognized by bowtie2).

First, be sure to move the Host_concat index genome to the same folder as the trimmed acer fasta files, and make sure they're not hidden in a subdirectory but directly in the same folder. 


```{bash}
#!/usr/bin/env bash
#BSUB -e bowtie2map_libro.err
#BSUB -o bowtie2map_libro.out
#BSUB -P and_transcriptomics
#BSUB -q general
#BSUB -n 8


cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer"

glob=".trim"

# Loop through files matching the pattern
for f in *$glob*; do
    output_file="${f%$glob}.sam"
    bowtie2 --local -f -U "$f" -x Host_concat --keep_unal -k 5 -S "$output_file"
done
```

Ok that code didn't work, "--keep_unal" isn't valid or recognized.

But in the error output file, the --un and --al flags do work, so include those instead.

```{bash}
#!/usr/bin/env bash
#BSUB -e bowtie2map_libro.err
#BSUB -o bowtie2map_libro.out
#BSUB -P and_transcriptomics
#BSUB -q general
#BSUB -n 8


cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer"

glob=".trim"

# Loop through files matching the pattern
for f in *$glob*; do
    output_file="${f%$glob}.sam"
    bowtie2 --local -f -U "$f" -x Host_concat --un --al -k 5 -S "$output_file"
done
```

Ok I got the same error as before: "Error: reads file does not look like a FASTA file
terminate called after throwing an instance of 'int'
(ERR): bowtie2-align died with signal 6 (ABRT) (core dumped)"

And I asked ChatGPT if one of my .trim files looked like a fasta file (by copying and pasting the first few lines in), and it DOESN'T look like a fasta file, but a fastq file. So that could be my issue. 

By default, bowtie2 takes a fastq file, so i just need to remove the "-f" flag (which tells it that it's fasta format).

Also, it wrote a file called "--al" so if i want the unaligned file to be written, i need to specify it's name/path.

Also, I don't know why I need the "--al" flag for writing unpaired reads that aligned at least once. Isn't that essentially what the sam file is? I'm gonna remove it.

```{bash}
#!/usr/bin/env bash
#BSUB -e bowtie2map_libro.err
#BSUB -o bowtie2map_libro.out
#BSUB -P and_transcriptomics
#BSUB -q general
#BSUB -n 8


cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer"

glob=".trim"

# Loop through files matching the pattern
for f in *$glob*; do
    output_file="${f%$glob}.sam"
    bowtie2 --local -U "$f" -x Host_concat --un "$output_file".unaligned -k 5 -S "$output_file"
done
```

Ok that's currently running so it might be working.

It ran and created .sam files and .sam.unaligned files for each Acer sample. Now I need to run the mapping efficiency code that Michael used and see how my numbers and his numbers matched up for Libro et al. 2013. 

```{perl}
#!/usr/bin/perl

print "

countreads_align.pl : counts the number of mapped reads in a bunch of fastq files
argument - glob to fastq files, default \.trim.al

";

my $glob="\.sam";
if ($ARGV[0]) { $glob=$ARGV[0];}

opendir THIS, ".";
my @fqs=grep /$glob/,readdir THIS;
my $f;
my $nrd;
foreach $f(@fqs){
	$nrd=`cat $f | wc -l`;
	$nrd=$nrd/4;
	print "$f\t$nrd\n";
}
```

countreads_align.pl > countreads_align.txt (ran this in the command line)

My numbers and Michael's numbers don't match at all. I am going to need him to show me what he did exactly.
