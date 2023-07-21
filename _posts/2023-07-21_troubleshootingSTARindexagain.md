---
layout: post
title: troubleshooting STAR index (again)
date: '2023-07-21'
categories: coding
tags: [coding]
---

Yesterday I had Natalia look at my code and we discovered something in the STAR index output (from the original [STAR index](https://github.com/ademerlis/ademerlis.github.io/blob/master/_posts/2023-06-21_STARindexAcergenome.md) and the [updated gff3 file STAR index after following Jill's code](https://github.com/ademerlis/ademerlis.github.io/blob/master/_posts/2023-06-27_stringtiecode_andredoingSTARwithupdatedannotationfile.md)). 

when I open the "geneInfo.tab" file from the STAR index, it looks like this:

<img width="305" alt="Screen Shot 2023-07-21 at 8 36 42 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/558e5bb5-4c3b-404d-9607-805e0d83a324">

and then when I look at the "Log.out" file from the STAR index, i see this warning: 

<img width="974" alt="Screen Shot 2023-07-21 at 8 45 25 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/e218e0e9-a331-4a5d-8999-c558ad3901b2">

So it seems like even with adding Jill's R script to edit the gff3 file, STAR still isn't able to recognize what is a gene id. Thus, all of them are "novel" and when I run stringtie, I get all novel gene names that stringtie makes up (Jill was the one who helped me notice this when I showed her my original counts matrix):

<img width="825" alt="Screen Shot 2023-07-21 at 8 47 56 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/a3c13e70-d89f-4788-847c-8aa2f8dd8167">

the "MSTRG" is a stringtie-generated gene ID, not one that is from any annotations from the gff3 or fasta file. Some of the acer ones mapped, but they are .pilon, which corresponds to exons. But the genes themselves do not have the proper IDs from the gff file. 

So what did we do? Well first we asked ChatGPT if "transcript_id" or "gene_id" are needed for STAR to index the genome properly. We think that is the step that failed because of the "geneInfo.tab" being blank. ChatGPT said:

"Yes, in a GFF3 (Generic Feature Format Version 3) file, the "gene_id" attribute is typically used to uniquely identify each gene feature. Including the "gene_id" attribute is essential, especially if you are working with gene annotations and performing downstream analyses that require tracking and identifying specific genes.

The "gene_id" attribute provides a unique identifier for each gene feature, which helps in linking together different features related to the same gene, such as exons, transcripts, and other regulatory elements. It simplifies data organization and facilitates querying and retrieval of information about specific genes during data analysis.

When creating or curating a GFF3 file, make sure that each gene feature (representing a gene locus) has a unique "gene_id" assigned to it. The value of "gene_id" can be any unique identifier, such as an alphanumeric code or a gene symbol, as long as it is unique across the entire GFF3 file.

In a GFF3 (Generic Feature Format Version 3) file, the "transcript_id" attribute is not strictly necessary for mRNA and exons, but it is highly recommended to include it for better data organization and downstream analysis.

The "transcript_id" attribute provides a unique identifier for each transcript feature. When included for mRNA and exons, it helps in linking these features together, indicating which exons belong to which transcript. This information is crucial for understanding the gene structure, alternative splicing events, and gene expression analysis.

Including the "transcript_id" attribute also makes it easier to track and retrieve information about specific transcripts during data analysis. It enables the association of exons with their corresponding mRNA transcripts, which is particularly important when studying gene isoforms and alternative splicing patterns.

While you can technically have mRNA and exons listed in a GFF3 file without the "transcript_id" attribute, doing so may lead to complications during downstream analyses, and it may be harder to interpret the data accurately. Therefore, it is considered good practice to include the "transcript_id" attribute for mRNA and exons in GFF3 files."

So essentially we need to manipulate the file in R so that the ID = gene_id so that STAR can properly index it.

So we tested this. First we ran this code in R, which basically adds on to Jill's coe and includes changing "ID=" to "gene_id=". 

```{r}
#Load libraries
library(tidyverse)

#Load  gene gff
Acerv.gff <- read.csv(file="Downloads/Galaxy1-[Acerv_assembly_v1.0.gff3].gff3", header=FALSE, sep="\t", skip=1) 

#rename columns
colnames(Acerv.gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
head(Acerv.gff)

# create gene id
Acerv.gff$gene_id <- gsub("ID=", "gene_id=", Acerv.gff$gene)

# Creating transcript id
Acerv.gff$transcript_id <- sub(";.*", "", Acerv.gff$gene)
Acerv.gff$transcript_id <- gsub("ID=", "", Acerv.gff$transcript_id) #remove ID= 

# delete old column
Acerv.gff %>% select(!gene) -> Acerv.gff 

# Checking what kinds of ids are in gff
unique(Acerv.gff$id)
# [1] "gene"        "mRNA"        "exon"        "CDS"         "start_codon" "stop_codon"  "tRNA"       

#If id == mRNA, exon, start_codon, stop_codon, CDS, tRNA, add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with original gene
Acerv.gff <- Acerv.gff %>% 
  mutate(gene_id = ifelse(id != "gene", paste0(gene_id, ";transcript_id=", Acerv.gff$transcript_id),  paste0(gene_id)))
head(Acerv.gff)

# Remove last col
Acerv.gff <- Acerv.gff[,-10]
head(Acerv.gff)  

#save file
write.table(Acerv.gff, file="~/Downloads/Acerv.GFFannotations.fixed_transcript.gff3", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)
```

Then, we tested whether this index worked by re-running the index code and then testing it on 4 samples. 

```{bash}
#!/bin/bash
#BSUB -J Acer_star_index_fixedannotations_take2
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -o /scratch/projects/and_transcriptomics/genomes/Acer/star_index%J.out
#BSUB -e /scratch/projects/and_transcriptomics/genomes/Acer/star_index%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

/scratch/projects/and_transcriptomics/programs/STAR-2.7.10b/bin/Linux_x86_64_static/STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeDir ${and}/genomes/Acer/Acer_STAR_index_gffannotations.fixed_take2 \
--genomeFastaFiles ${and}/genomes/Acer/Acerv_assembly_v1.0_171209.fasta \
--sjdbGTFfile ${and}/genomes/Acer/Acerv.GFFannotations.fixed_transcript_take2.gff3 \
--sjdbOverhang 100 \
--sjdbGTFtagExonParentTranscript Parent \
--genomeSAindexNbases 13
```
And the resulting geneInfo.tab looks like this:

<img width="580" alt="Screen Shot 2023-07-21 at 8 58 11 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/7274a0b5-c4b7-43c0-a079-4aa8aec20658">

I asked ChatGPT whether "MissingGeneType" was important or not, and it said:

"In the context of the STAR aligner, "MissingGeneType" is a warning message that might appear during the genome indexing process. This warning is not an error that would prevent STAR from successfully indexing the genome, but it indicates that some gene types were not found or recognized during the indexing step.

The warning message "MissingGeneType" occurs when the gene type information is provided in the GTF (Gene Transfer Format) annotation file used for genome indexing, but certain gene types are not specified or not recognized by STAR. The gene types in the GTF file are typically specified in the third field of each GTF entry (columns are tab-separated), which is the "feature" field.

Here's an example of a GTF entry with the "feature" field representing the gene type:

```
chr1  Ensembl gene  11869  14409  .  +  .  gene_id "ENSG00000223972"; gene_name "DDX11L1";
```

If STAR encounters gene types in the GTF that it does not recognize (e.g., "mRNA", "exon", etc.), it will display the "MissingGeneType" warning to alert you that certain gene types were not found in the pre-defined list of supported gene types.

While this warning message should not affect the indexing process itself, it's essential to review your GTF file and ensure that the gene types you use are compatible with STAR's expected gene type names. Additionally, make sure that the GTF file is correctly formatted and follows the GTF specification. Correcting any issues in the GTF file will help ensure the accuracy of downstream analysis when using STAR for read alignment against the indexed genome."


Ok, when I was looking through my tabs, I saw Ariana Huffmyer's code had this file: "fix_gff_format.Rmd". When I look at it, she has the code in there to add both the transcript_id and the gene_id to the gff file. While she is using a different species, it has the same format as the Acer gff3 file I'm using ("ID=" and "Parent=").

Her code does this:

```{r}
---
title: "Fix GFF format"
author: "Ariana S Huffmyer"
date: "2022"
output: html_document
editor_options: 
  chunk_output_type: console
---
#This script add transcript and gene id into GFF file for alignment.  

Here, I'll be adding transcript_id= and gene_id= to 'gene' column because we needs that label to map our TagSeq data  

Load libraries and data. 

#Load libraries
library(tidyverse)
library(R.utils)


Load  gene gff file

gff <- read.csv(file="Mcap2020/Data/TagSeq/Montipora_capitata_HIv3.genes.gff3", header=FALSE, sep="\t") 

Rename columns 

colnames(gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
head(gff)

Create transcript ID  

gff$transcript_id <- sub(";.*", "", gff$gene)
gff$transcript_id <- gsub("ID=", "", gff$transcript_id) #remove ID= 
gff$transcript_id <- gsub("Parent=", "", gff$transcript_id) #remove ID= 
head(gff)


Create Parent ID 

gff$parent_id <- sub(".*Parent=", "", gff$gene)
gff$parent_id <- sub(";.*", "", gff$parent_id)
gff$parent_id <- gsub("ID=", "", gff$parent_id) #remove ID= 
head(gff)


Now add these values into the gene column separated by semicolons.  

gff <- gff %>% 
  mutate(gene = ifelse(id != "gene", paste0(gene, ";transcript_id=", gff$transcript_id, ";gene_id=", gff$parent_id),  paste0(gene)))
head(gff)


Now remove the transcript and parent ID separate columns.  


gff<-gff %>%
  select(!transcript_id)%>%
  select(!parent_id)

head(gff)


Save file. Then upload this to Andromeda for use in bioinformatic steps.  

write.table(gff, file="Mcap2020/Data/TagSeq/Montipora_capitata_HIv3.genes_fixed.gff3", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

#gzip the file 
gzip("Mcap2020/Data/TagSeq/Montipora_capitata_HIv3.genes_fixed.gff3")
```

When I look at my most updated gff3 file ("Acerv.GFFannotations.fixed_transcript_take2.gff3"), I see that "Parent=" is still there. I don't know if this messes it up at all, since I think "Parent" is supposed to correspond to the gene_id, and then the "transcript_id" is supposed to be the more specific one (i.e. the one that has the .exon1 in it). 

<img width="1266" alt="Screen Shot 2023-07-21 at 9 33 17 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/dd5a48e7-7152-438a-880f-30443e8ea4e9">

Just to be safe, I'm going to try indexing with a 3rd version of the gff3 file following Ariana's Rmd code. Because it looks weird to me that my "gene_id" ones in the gff3 file have the ".exon" in their names. I don't think that's right.

Here's my most updated gff fixing R code:

```{r}
#Load libraries
library(tidyverse)

#Load  gene gff
Acerv.gff <- read.csv(file="Downloads/Galaxy1-[Acerv_assembly_v1.0.gff3].gff3", header=FALSE, sep="\t", skip=1) 

#rename columns
colnames(Acerv.gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
head(Acerv.gff)

# Creating transcript id
Acerv.gff$transcript_id <- sub(";.*", "", Acerv.gff$gene)
Acerv.gff$transcript_id <- gsub("ID=", "", Acerv.gff$transcript_id) #remove ID= 
Acerv.gff$transcript_id <- gsub("Parent=", "", Acerv.gff$transcript_id) #remove Parent=
head(Acerv.gff)

# Create Parent ID 
Acerv.gff$parent_id <- sub(".*Parent=", "", Acerv.gff$gene)
Acerv.gff$parent_id <- sub(";.*", "", Acerv.gff$parent_id)
Acerv.gff$parent_id <- gsub("ID=", "", Acerv.gff$parent_id) #remove ID= 
head(Acerv.gff)

Acerv.gff <- Acerv.gff %>% 
  mutate(gene = ifelse(id != "gene", paste0(gene, ";transcript_id=", Acerv.gff$transcript_id, ";gene_id=", Acerv.gff$parent_id),  paste0(gene)))
head(Acerv.gff)

Acerv.gff<-Acerv.gff %>%
  select(!transcript_id)%>%
  select(!parent_id)

#save file
write.table(Acerv.gff, file="~/Downloads/Acerv.GFFannotations.fixed_transcript_take3.gff3", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)
```

Now it looks like this. This might be right? Idk how much this changes things.

<img width="1337" alt="Screen Shot 2023-07-21 at 9 46 43 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/13946726-4337-41c4-ae12-e384caff488a">

Let's try now running the alignment on the samples. NOTE: I have two pipelines of trimmed sequences now: the NJ_pipeline is from Natalia and Jill and that involved using cutadapt. The AS_pipeline is Ariana and Sam and involved using fastp. I am going to try this new STAR gff file on the NJ pipeline first since that is farther along (all the way to counts matrix) and so it can be directly compared.

