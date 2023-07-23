---
layout: post
title: updated stringtie code
date: '2023-07-22'
categories: coding
tags: [coding, Ch4_AcerCCC]
---

Following [Sam's code](https://github.com/SamGurr/SamGurr.github.io/blob/master/_posts/2021-01-07-Geoduck-TagSeq-Pipeline.md#Trimming-polyA-tail), I used these arguments below for stringtie:

Main StringTie arguments used below:

-p = specify the number of threads (CPUs). Note that a single node has at least 8 CPUs or 'threads' to call here

-A = Gene abundances will be reported (tab delimited format) in the output file with the given name.

-e = Limits the processing of read alignments to only estimate and output the assembled transcripts matching the reference transcripts given with the -G option (requires -G, recommended for -B/-b). With this option, read bundles with no reference transcripts will be entirely skipped, which may provide a considerable speed boost when the given set of reference transcripts is limited to a set of target genes, for example.

 -G <ref_ann.gff> = Use the reference annotation file (in GTF or GFF3 format) to guide the assembly process. The output will include expressed reference transcripts as well as any novel transcripts that are assembled. This option is required by options -B, -b, -e, -C.

-o = Sets the name of the output GTF file where StringTie will write the assembled transcripts

```{bash}
#!/bin/bash
#BSUB -J stringtie_updatedannotations_take3
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o stringtie_updatedannotations_take3%J.out
#BSUB -e stringtie_updatedannotations_take3%J.err
#BSUB -u and128@miami.edu
#BSUB -N

module load python/3.8.7

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Ch4_AcerCCC/aligned_updatedannotations_take3"

data=($(ls *Aligned.sortedByCoord.out.bam))

for i in ${data[@]} ;

do \
/scratch/projects/and_transcriptomics/programs/stringtie-2.2.1/stringtie -p 8 -e -B -G /scratch/projects/and_transcriptomics/genomes/Acer/Acerv.GFFannotations.fixed_transcript_take3.gff3 -A ${i}.gene_abund.tab -o ${i}.gtf ${i} ; \
done
```


