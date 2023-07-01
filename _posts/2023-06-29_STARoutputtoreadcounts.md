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
