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

36455 sequences.
1678 average length.
65308 maximum length.
61 minimum length.
N50 = 2747
61.2 Mb altogether (61163689 bp).
0 ambiguous Mb. (100 bp, 0%)
0 Mb of Ns. (100 bp, 0%)


Next step is to get the uniprot annotations with blast.

```{bash}
# getting uniprot_swissprot KB database
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

gunzip uniprot_sprot.fasta.gz
```

Then I created a script for the makedatabase part.

```{bash}
#!/usr/bin/bash
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -W 120:00
#BSUB -o makedatabase.out
#BSUB -e makedatabase.err

module load blast/2.2.29+
makeblastdb -in uniprot_sprot.fasta -dbtype prot
```


That ran. 

Now, I need to split up the fasta file so the blast code runs quickly.

```{bash}
#navigate to wd where splitFasta.pl is located 
#copy fasta file into this directory
cp /scratch/projects/and_transcriptomics/genomes/Acer/Locatelli_2023/Acer_Genome/Acropora_cervicornis.mrna-transcripts.fa Acer2023.fasta

#run split Fasta code
perl splitFasta.pl Acer2023.fasta 200
```

This created 200 subset.fasta files. Now, need to create a script to submit a job to pegasus that blasts all 200 chunks to uniprot in parallel.

I think I'll need to create a for-loop so it creates individual jobs for each chunk, so that they can run simultaneously. 

This for loop worked: 

```{bash}
#BSUB -u and128@miami.edu

#specify variables and paths

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="/scratch/projects/and_transcriptomics/genomes/Acer/annotatingTranscriptomes-master"

cd ${projdir}

data=($(ls subset*))

for file in "${data[@]}" ; do \

#build script
echo "making blast uniprot script for ${file}..."
echo "
#! /usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${file}_blast_uniprot
#BSUB -e ${projdir}/logs/${file}_blast_uniprot.err
#BSUB -o ${projdir}/logs/${file}_blast_uniprot.out
#BSUB -W 12:00
#BSUB -n 8
#BSUB -q general

cd ${projdir}

module load blast/2.2.29+

blastx -query ${file} -db uniprot_sprot.fasta -evalue 0.0001 -num_threads 4 -num_descriptions 5 -num_alignments 5 -out ${file}.br" > ${projdir}/${file}_blast_uniprot.job

bsub < ${projdir}/${file}_blast_uniprot.job

done
```

This took ~ l hr to complete.

To make sure it worked, run the following line and you should end up with the same number of queries as there are sequences from the seqstats_Acer.txt file (created from running seq_stats.pl). 

```{bash}
grep "Query= " subset*.br | wc -l
```

I get 36,455, which equals the number of sequences in the seqstats_Acer.txt file.


Next, combine all blast results:

```{bash}
cat subset*br > myblast.br
```

Then, extract coding sequences and corresponding protein translations using the perl script **CDS_extractor_v2.pl**. 

This is a really long script so I'm not going to paste it in here. I'm first going to try to run it as-is. It looks like it requires Bio::SeqIO so we'll see how that goes.

I should probably create a job script and submit this as a job just in case it takes awhile.

Michael runs CDS_extractor_v2.pl on a file he made called "Acervicornis_iso.fasta". It looks like he ran this code on the original fasta file:

```{bash}
cat Acervicornis.fasta | perl -pe 's/>Acropora(\d+)(\S+).+/>Acropora$1$2 gene=Acropora$1/'>Acervicornis_iso.fasta
```

He had flagged this line of code for only being specific to Trinity-assembled transcriptomes.

The one I'm using was generating using "funannotate." So, each header of the fasta file looks something like this:

```{bash}
>FUN_000006-T1 FUN_000006
```

I can try to adapt Michael's code so it does the same thing, which is basically substituting the trinity header to be named an arbitrary gene name (i.e. Acropora00001). 

**update**: actually i'm not sure I want to do that, because then I would also have to edit the enitre gff3 file (the gene IDs in there are FUN_#). I'm going to leave it as is and change the Acervicornis_seq2iso.tab code to grep and awk FUN instead of Acropora.

Remember the Acer_2023.fasta is the mrna_transcripts.fa file that I just copied over from the original genome download folder.

"Acervicornis_seq2iso.tab"  is the necessary file for generating read counts per gene back in the tagSeq bioinformatics pipeline (which is the whole reason i'm running this annotatingTranscriptomes pipeline).

So i should run both lines to get both files: Acervicornis_seq2iso.tab and Acervicornis_iso.fasta.

```{bash}
grep ">" Acer2023.fasta | perl -pe 's/>FUN(\d+)(\S+)\s.+/FUN$1$2\tFUN$1/'>Acervicornis_seq2iso.tab

cat Acer2023.fasta | perl -pe 's/>FUN(\d+)(\S+).+/>FUN$1$2 gene=FUN$1/'>Acervicornis_iso.fasta
```

Ok the second line doesn't do anything because the original fasta file is already listed in this exact format. So i'm not sure if any of this code is actually useful.


```{bash}
# extracting coding sequences and corresponding protein translations:
conda activate bioperl #if not active already
echo "perl ~/bin/CDS_extractor_v2.pl Acervicornis_iso.fasta myblast.br allhits bridgegaps" >cds

#! /usr/bin/env bash
#BSUB -P and_transcriptomics
#BSUB -J CDS_extractor
#BSUB -e /scratch/projects/and_transcriptomics/genomes/Acer/annotatingTranscriptomes-master/CDS_extractor.err
#BSUB -o /scratch/projects/and_transcriptomics/genomes/Acer/annotatingTranscriptomes-master/CDS_extractor.out
#BSUB -W 12:00
#BSUB -n 8
#BSUB -q general

cd "/scratch/projects/and_transcriptomics/genomes/Acer/annotatingTranscriptomes-master"

perl CDS_extractor_v2.pl Allyson_Created_Files/Acer2023.fasta Allyson_Created_Files/myblast.br allhits bridgegaps

```

To run the above job script, I had to manually install the perl module "Bio::DB::Fasta:. I ran this:

```{bash}
cpan Bio::DB::Fasta
```
But this doesn't work because there are a ton of other dependencies it needs. I'm going to abandon this pipeline at this step because there is an annotation file that came with the Locatelli 2023 Acer genome and it has GO terms in it. I'll just move forward with that for right now.

