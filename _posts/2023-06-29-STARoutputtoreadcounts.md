---
layout: post
title: STAR output to gene counts
date: '2023-06-29'
categories: coding
tags: [coding, CCC_ch4]
---

I can't get stringtie to work so I want to try Natalia's method instead, where she took the STAR read counts and somehow turned that into a gene count matrix. 

First, she downloads the multiqc_general_stats.txt file from running multiqc on the aligned STAR reads. Then she does a lot of data tidying in R ([code here](https://github.com/China2302/SCTLD_RRC/blob/main/05_read_counts.Rmd)

```{r}
library(data.table)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(reshape2)
library(tidyverse)
```


```{r}
#Sorting samples with less than 4M and 5M reads mapped

reads_mapped<- read.delim("multiqc_general_stats.txt", sep = "\t")
head(reads_mapped)

reads_mapped<- reads_mapped %>% 
           rename(percent_uniq = STAR_mqc.generalstats.star.uniquely_mapped_percent) %>%
  
               rename(million_uniq = STAR_mqc.generalstats.star.uniquely_mapped) %>% 
  
               rename(Sample_ID = Sample)

samples_selected_4M<- reads_mapped %>% filter(million_uniq > 4000000)

samples_selected_5M<- reads_mapped %>% filter(million_uniq > 5000000)
```


```{r}
# read in metadata
sample_metadata <- read.csv("sample_metadata.csv")
names(sample_metadata)

#join read counts
sample_metadata_reads<- sample_metadata %>% left_join(samples_selected_4M, by="Sample_ID") 

sample_metadata_reads <- drop_na(sample_metadata_reads)
```


```{r}
# create a list of all files from samples that have at least 4M reads mapped

file_4M<- list.files("star_trimmed_reads",
                  "*.gzReadsPerGene.out.tab", full.names = T)

countData_4M = data.frame(fread(file_4M[1]))[c(1,3)]

for(i in 2:length(file_4M)) {
        countData_4M = cbind(countData_4M, data.frame(fread(file_4M[i]))[3])
}

# Skip first 4 lines, count data starts on the 5th line
countData_4M = countData_4M[c(5:nrow(countData_4M)),]
colnames(countData_4M) = c("GeneID", gsub(paste0(dir,"samples_4M_reads/"), "", file_4M))
colnames(countData_4M) = gsub("_trimmed.fq.gzReadsPerGene.out.tab", "", colnames(countData_4M))
rownames(countData_4M) = countData_4M$GeneID

countData_4M = countData_4M[,c(2:ncol(countData_4M))]

#write_rds(countData_4M, "data/countData_4M.rds")
```

Ok, when I get to this part of the code, I realize that my files look different than Natalia's. Her column 5 in the countData_4M has gene counts, while mine says "missing gene ID". 

<img width="192" alt="Screen Shot 2023-07-01 at 3 41 34 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/14be7712-1139-40a1-8e5f-4e8278fb7c72">


So when I follow all her code in the last code chunk, I end up with this for all samples:

<img width="906" alt="Screen Shot 2023-07-01 at 3 42 11 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/978ce27e-fcce-4e8c-8bee-474db977e5b1">

But in her code it looks like she is able to name the rownames as "Gene ID", meaning that the rows correspond to genes and the values under each sample are the counts. 

So what does this mean? did she use STAR to do the quantreads step? Or does this mean my alignmnet didn't work with the annotations? 

I just checked several files for both the aligned and anilded_updatedannotations files (*ReadsperGene.out.tab) and they all have the "missing gene id" thing. Then I looked at the other out.tab file, the one with "SJ" which means splice junctions, and that had a ton of genes listed there. So maybe nothing aligned to the reference genome and actually everything is being considered as "novel" which is why they are read as splice junctions (i think?). 

I asked ChatGPT what it could mean, and this is what it said: 

"If you are seeing "missing gene ID" entries in the `readsPerGene.tab` file generated by the STAR aligner for RNA-seq data, it typically indicates that the reads from your sequencing data did not align to any known gene in the reference genome used for alignment.

Here are a few possible reasons for this:

1. Mismatched reference genome: The reference genome used for alignment might not match the organism or the specific version of the genome you are working with. It is crucial to ensure that you are using the correct reference genome for your analysis.

2. Novel genes or isoforms: If the reference genome you used does not contain all the genes or isoforms present in your specific sample or experimental condition, it is possible that some reads align to regions that are not annotated in the reference genome. These reads will be labeled as "missing gene ID" in the `readsPerGene.tab` file.

3. Low-quality reads: It's possible that the reads with "missing gene ID" did not align to any genes due to low-quality sequence data, sequencing errors, or other technical issues. In such cases, it is recommended to check the quality of your reads, perform any necessary preprocessing steps, and ensure that your sequencing data is of sufficient quality.

To further investigate the issue, you can examine other alignment-related files generated by STAR, such as the alignment summary, log files, or BAM files, which may provide additional information about the alignment process and the mapping results.

Additionally, if you suspect the absence of certain genes or isoforms in the reference genome, you might consider using alternative alignment tools, such as HISAT2 or Salmon, which are specifically designed to handle novel transcripts or splice junctions that are not present in the reference genome."

Ok, my first thought is maybe that the samples/sequences are too low-quality and that's why they aren't aligning properly. But I guess waht Im stuck on is if they didn't align properly, then why do the STAR outputs give a percentage of mapped reads? to me that means the alignmnet worked, but the gene ids don't match up for some reason. Does that mean there is something wrong with the genome index? the .gff3 file has the gene names in it, so if a read aligned to it why wouldnt it be annotated. 

This is what ChatGPT says:
"If you are getting a percentage of reads mapped in STAR but also encountering "missing gene IDs" in the results, it suggests that while a proportion of your reads have successfully aligned to the reference genome, they might not be aligning to any annotated genes present in the reference genome.

Here are a few potential reasons for this situation:

Unannotated or novel genes: It's possible that some of the aligned reads correspond to unannotated or novel genes that are not included in the reference genome. These genes might not have official gene IDs or annotations, resulting in "missing gene IDs" in the results.

Non-coding RNA or regulatory elements: RNA-seq reads might align to non-coding RNA molecules or regulatory elements, such as long non-coding RNAs or enhancers, which might not have gene IDs associated with them in the reference genome annotation.

Mapping to repetitive regions or non-genic regions: Some reads may align to repetitive regions of the genome or regions that do not correspond to annotated genes. These regions might not have gene IDs associated with them, resulting in "missing gene IDs" for those reads.

Low-quality alignments: While a percentage of reads may be mapped, it's possible that some of these alignments are of lower quality. They might align to incorrect regions or have a higher mapping ambiguity, leading to incorrect or missing gene annotations.

To further investigate this issue, you can examine the alignment files generated by STAR, such as the BAM file, and explore the specific alignments for the reads with "missing gene IDs". This can help identify the regions or features to which those reads are aligning and provide insights into why they might not have associated gene IDs.

Additionally, considering alternative alignment tools or exploring specialized pipelines for identifying novel genes, non-coding RNAs, or regulatory elements could provide further insights into your data."
