---
published: false
---
## Finally getting the gene counts for Acer samples from Chapter 2

Following the last section of [Michael Studivan's code](https://github.com/mstudiva/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt). 

```{bash}
# NOTE: Must have a tab-delimited file giving correspondence between contigs in the transcriptome fasta file and genes

cat Acer/Locatelli_2023/Acervicornis_seq2iso.tab Symbiodinium/Mstudiva_Symbiodinium_annotated_transcriptome_code/Symbiodinium_seq2iso.tab > Host_concat_seq2iso.tab
```

Job: use sam tools and [samcount.pl](https://github.com/mstudiva/tag-based_RNAseq/blob/master/samcount.pl).

arg1: SAM file (by cluster, contig, or isotig)
arg2: a table in the form 'reference_seq<tab>gene_ID', giving the correspondence of 
reference sequences to genes. With 454-deived transcriptome, the gene_ID would be isogroup; 
with Trinity-derived transcriptiome,it would be component.
  
paths:
  1. /scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer/bowtie2align_LocatelliShoguchi/sam_files/*.sam
  2. /scratch/projects/and_transcriptomics/genomes/Host_concat_seq2iso.tab
  
I need to write a for-loop to each .sam file to do this: print "samcount.pl $file $tab aligner=bowtie2 >$file.counts" but written in bash not perl.
  
```{bash}
  
```
  
  
  








