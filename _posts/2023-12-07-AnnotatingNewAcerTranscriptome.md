---
layout: post
title: Annotating new Acer transcriptome
date: '2023-12-07'
categories: Coding
tags:
  - Coding
published: true
---

I'm going to be following [Michael Studivan's code](https://github.com/mstudiva/Acropora-cervicornis-annotated-transcriptome/blob/main/tagSeq_TranscriptomeAnnotation_README.txt) for annotating the new Acer transcriptome. 

Step 1 is installing bioperl. I had trouble doing this, but I didn't realize I could activate it as a conda environment. It looks like it's working so far.

```{bash}
# To install Bioperl as a conda environment
conda create -y -n bioperl perl-bioperl

conda activate bioperl

git clone https://github.com/z0on/annotatingTranscriptomes.git

git clone https://github.com/z0on/emapper_to_GOMWU_KOGMWU.git

seq_stats.pl /scratch/projects/and_transcriptomics/genomes/Acer/Locatelli_2023/Acer_Genome/Acropora_cervicornis.mrna-transcripts.fa  > seqstats_Acer.txt

```

I had several issues running the seq_stats.pl script.

First, I had to update the shebang in the .pl script so that it was calling the right perl (within the bioperl environment).

Then, I also had to remove the line `use Bio::Perl;` because it was causing issues.

So the seq_stats.pl script beginning looks like this now:

```{perl}
#!/usr/anaconda3/envs/bioperl/bin/perl

use lib "~/anaconda3/envs/bioperl";
use Bio::SeqIO;
```

Now I run this with the updated script:

```{bash}
seq_stats.pl /scratch/projects/and_transcriptomics/genomes/Acer/Locatelli_2023/Acer_Genome/Acropora_cervicornis.mrna-transcripts.fa  > seqstats_Acer.txt
```
Results: 
-------------------------
36455 sequences.
1678 average length.
65308 maximum length.
61 minimum length.
N50 = 2747
61.2 Mb altogether (61163689 bp).
0 ambiguous Mb. (100 bp, 0%)
0 Mb of Ns. (100 bp, 0%)
-------------------------


Next step is to get the uniprot annotations with blast.

```{bash}
# getting uniprot_swissprot KB database
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
```

