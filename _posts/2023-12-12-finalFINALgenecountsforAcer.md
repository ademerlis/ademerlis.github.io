---
published: false
---
## Finally getting the gene counts for Acer samples from Chapter 2

Following the last section of [Michael Studivan's code](https://github.com/mstudiva/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt). 

```{bash}
# NOTE: Must have a tab-delimited file giving correspondence between contigs in the transcriptome fasta file and genes

cat Acer/Locatelli_2023/Acervicornis_seq2iso.tab Symbiodinium/Mstudiva_Symbiodinium_annotated_transcriptome_code/Symbiodinium_seq2iso.tab > Host_concat_seq2iso.tab





