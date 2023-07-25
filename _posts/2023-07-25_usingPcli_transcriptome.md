---
layout: post
title: using Pcli transcriptome
date: '2023-07-25'
categories: coding
tags: [coding, Ch2_tempvariability]
---

So while the stringtie assembly is running for Acer, I am going to try to figure out how to use the Pcli file I downloaded (see [previous blog post](https://github.com/ademerlis/ademerlis.github.io/blob/master/_posts/2023-07-21_Pclitranscriptome.md)). 

First i gunzipped it so it would be .fna format. 

Then I just did "head" to see what the file looked like:

<img width="690" alt="Screen Shot 2023-07-25 at 9 59 24 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/36eafd47-8ef3-42b8-942b-337f2220eb55">

the "len=" and "path=[...]" are not normally part of a fasta file. 

So I'm confused about many things related to this -- so I have this transcriptome that was assembled de novo from a metatranscriptome by Avila-Magana et al. 2021. But we just have the .fna file, no gff file or anything. In reading their methods, it sounds like they used Kallisto to do pseudoalignment to a de novo reference transcriptome (so using their own sequences and then tring to align to that?) and then they used GO_MWU to get gene ontology terms. But how did they go from the transcriptome to GO_MWU, because that requires gene IDs with GO terms. 

Ok in their methods they say this: "Amino acid sequences were predicted from the Coral host and Symbiodiniaceae transcriptomes by using Transdecoder (v2.0.1)87. Orthologous groups of protein sequences amongst the three coral and their associated three photosymbiont species were determined with the OrthoFinder (v2.4.0)16 bioinformatics tool, using default parameters. Using reciprocal best-hits via BLAST all-v-all algorithm, Orthofinder determined the number of conserved putative orthologues among the three coral and three photosymbiont species. By a custom python script, the ortholog differential gene expression (logFC) was retrieved for each expression dataset per species against the other two species."

That sounds like how they did their "annotations". Because what I'm missing is basically the connection from the host transcriptome to any sort of identifiable information. 

Or what is we do something like this:

" To generate gene ontology annotations from our 169 assemblies we utilized the annotation software EggNog mapper (Huerta-Cepas et al., 2017, 170 2019) through the online web portal against the eukaryote database. To annotate individual 171 contigs we used blastp against the swiss-prot uniprot database (The UniProt Consortium et al., 172 2021)." 
(from [Dimos et al. 2022](https://doi.org/10.1101/2021.04.28.441826)).

Is that how you make a gff file? 

Let's try using eggNog mapper to annotate the transcriptome, because that's also what it says they did in Avila-Magana et al. 2021 (before the GO-MWU step).

WELP. after all of this researching, I found the emapper annotation file on the dryad repository for the Avila-Magana et al. 2021 publication. https://datadryad.org/stash/dataset/doi:10.5061/dryad.k3j9kd57b

So this was all a huge waste of time. oh well.

So now to do step 1 of the Pcli stuff, which is the HISAT2 code.

First off, we have to edit the Dip_Host.fna file (which is a fasta file) because when it was assembled in Trinity, it added trinity-specific things in there that we don't need (and HISAT2 will get tripped up by).

<img width="798" alt="Screen Shot 2023-07-25 at 1 09 11 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/6f62367f-310c-4958-8ede-fa78cc0d4977">

According to ChatGPT, the "len=" and "path=" are:

"In the context of Trinity, the .fna file contains the assembled transcript sequences in FASTA format, along with additional annotations specific to Trinity's assembly output. The annotation information you provided appears to be related to the contig connectivity and pathing in the de Bruijn graph-based assembly.

Here's what the annotations in the .fna file likely mean:

len=261: This represents the length of the assembled transcript sequence, which is 261 nucleotides long.

path=[812:0-112 813:113-121 814:122-136 815:137-145 816:146-157 817:158-186 818:187-260]: This part corresponds to the pathing information of the contig in the de Bruijn graph. It indicates the connectivity of the contig in the graph and the positions of different segments within the contig. Each segment is defined by the node IDs in the graph (e.g., 812, 813, etc.) and their respective start and end positions within the contig (e.g., 0-112, 113-121, etc.).

[-1, 812, 813, 814, 815, 816, 817, 818, -2]: The numerical values here are node IDs representing the start and end of the contig in the de Bruijn graph. The values -1 and -2 indicate the start and end nodes, respectively.

Trinity's assembly output includes these pathing and connectivity annotations to help describe how the contigs were constructed from the raw reads during the de novo transcriptome assembly process.

If you are using this Trinity-generated .fna file as a reference for read alignment or downstream analysis, you typically don't need to parse or interpret these annotations directly. Instead, you can use the assembled transcript sequences as-is for read mapping, differential expression analysis, or other transcriptomic analyses. Trinity usually provides the primary assembly file (e.g., Trinity.fasta) that contains the sequences without the additional pathing information, and that file is often used for further analyses."

Of course though on the Dryad repository they only include the .fna file, not the Trinity.fasta. So I'll need to manually remove these things.

Thankfully ChatGPT also had a suggestion for that:

"To use Trinity's transcriptome assembly output with HISAT2, you will need to generate a clean version of the reference transcriptome in standard FASTA format. You can do this by removing the len= and path= annotations from the original .fna file. Typically, Trinity provides a primary assembly file, such as Trinity.fasta, which contains the assembled transcript sequences without the additional annotations. This primary assembly file is the one you should use for downstream analyses, including read alignment with HISAT2.

To convert the Trinity .fna file to a clean FASTA format suitable for HISAT2, you can use standard text-processing tools or scripting languages to parse the file and extract only the sequence information while discarding the len= and path= lines. Here's a simple example of how you can achieve this using the awk command in Linux:
```{bash}
# Assuming the Trinity-generated .fna file is named "transcriptome.fna"
awk '/^>/ {print $1}' transcriptome.fna > clean_transcriptome.fasta
```

So i'll make a new Dip_Host.fna file. 

Ok all that did was extract the gene/transcript IDs... and got rid of the sequences.

Ok try this instead:
```{bash}
awk '/^>/ {print; getline; print}' transcriptome.fna > clean_transcriptome.fasta
```

Ok that didn't work either. I had to go back and forth with ChatGPT a few times before we got one that worked.

```{bash}
sed '/^>/ s/ len=[0-9]* path=\[[^]]*\]//' transcriptome.fna > clean_transcriptome.fasta
```

This got us much closer. I just need to remove the extra string now that comes after the path=[...]. it also is [...]

```{bash}
sed '/^>/ s/ \[.*\]//' clean_Pcli_transcriptome.fasta > clean_Pcli_transcriptome_final.fasta
```

Alrighty it worked!!

<img width="371" alt="Screen Shot 2023-07-25 at 1 24 50 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/05b2c95b-b098-4ec3-aa9d-29f2a3460145">

So the updated code for Pcli HISAT2 will be:

```{bash}
#!/bin/bash
#BSUB -J HISAT2_Pcli
#BSUB -q bigmem
#BSUB -n 16
#BSUB -P and_transcriptomics
#BSUB -o HISAT2_Pcli%J.out
#BSUB -e HISAT2_Pcli%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

module load samtools/1.3
module load python/3.8.7

/scratch/projects/and_transcriptomics/programs/hisat2-2.2.1/hisat2-build -f ${and}/genomes/Pcli/clean_Pcli_transcriptome_final.fasta ${and}/genomes/Pcli/Pcli_reference_genome_hisat2
echo "Reference genome indexed. Starting alignment" $(date)

cd /scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/AS_pipeline/3_trimmed_fastq_files/Pcli_fastq_files/
array=($(ls *.fastq.gz))
for i in ${array[@]};
 do \
     	sample_name=`echo $i| awk -F [.] '{print $2}'`
        /scratch/projects/and_transcriptomics/programs/hisat2-2.2.1/hisat2 -p 8 --dta -x ${and}/genomes/Pcli/Pcli_reference_genome_$
        samtools sort -@ 8 -o ${sample_name}.bam ${sample_name}.sam
                echo "${i} bam-ified!"
        rm ${sample_name}.sam ;
done
```

But do I need to make any changes to the HISAT2 code since this is a transcriptome, not a genome?

Ok I asked ChatGPT:


To specify that you are using a transcriptome instead of a genome with HISAT2, you need to provide a transcriptome index during the indexing step and use the appropriate options during the alignment step.

Here's how you can create a transcriptome index and perform alignment with HISAT2:

Indexing the Transcriptome:
To create a HISAT2 transcriptome index, you need to use the hisat2-build command with your transcriptome reference in FASTA format. The command would look like this:

bash
Copy code
# Assuming your transcriptome reference is in "transcriptome.fasta" and you want to create an index named "transcriptome_index"
hisat2-build transcriptome.fasta transcriptome_index
This command will generate a set of files with the prefix "transcriptome_index" that constitute the HISAT2 index for your transcriptome.

Alignment with HISAT2:
To align your RNA-seq reads to the transcriptome using the HISAT2 index, you should use the hisat2 command with the appropriate options. Here's a basic example:

```{bash}
# Assuming your HISAT2 index is named "transcriptome_index" and your paired-end reads are in "read1.fastq" and "read2.fastq"
hisat2 -x transcriptome_index -1 read1.fastq -2 read2.fastq -S alignment.sam
```

In this example:
-x transcriptome_index: Specifies the HISAT2 index to be used for alignment, which is the transcriptome index you created earlier.
-1 read1.fastq -2 read2.fastq: Specifies the paired-end read files to be aligned to the transcriptome.
-S alignment.sam: Specifies the output file in SAM format, which will contain the alignment results.


**Please note that HISAT2 is primarily designed for genome alignment, and using it with a transcriptome has some limitations. Since transcriptomes do not contain introns, you may not get accurate information about splicing events during alignment. For better analysis of RNA-seq data, especially if you are interested in transcript isoforms and splicing, specialized tools such as STAR or kallisto that are designed for transcriptome alignment and quantification may be more suitable choices.**"

Ok so for transcriptome alignment, I shouldn't be using HISAT2. Maybe this is when I can try to use salmon?
