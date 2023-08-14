---
layout: post
title: Ch2 Acer DESeq2
date: '2023-08-14'
categories: coding
tags: [coding, ch2_tempvariability]
---

So I started on the code to get through DESeq2 for the ch2 temperature variability samples. 

First, I needed to tidy up the metadata file and figure out how many samples had at least 4 million reads post-alignment step.

## 1. obtaining_aligned_reads.Rmd

Graphs of alignment rates and reads aligned
```{r}
sequencing_data <- read.csv("../../RNA_extraction_sequencing_data.csv")

sequencing_data %>% 
  select(Species:Treatment,Raw.Reads..Million.,M.Reads.After.Filtering..Trimmed.reads.,Mapping.Rate..Post.trimming.alignment.) %>% 
  mutate(percent_alignment = gsub("%","", Mapping.Rate..Post.trimming.alignment.)) %>% 
  mutate(percent_alignment = as.numeric(percent_alignment)) %>% 
  ggplot(., aes(x=Sample.ID, y=percent_alignment, fill=Treatment)) + geom_col() + facet_wrap(~Species) + theme_classic()

sequencing_data %>% 
  select(Species:Treatment,Raw.Reads..Million.,M.Reads.After.Filtering..Trimmed.reads.,Mapping.Rate..Post.trimming.alignment.) %>% 
  mutate(percent_alignment = gsub("%","", Mapping.Rate..Post.trimming.alignment.)) %>% 
  mutate(percent_alignment = as.numeric(percent_alignment)) %>% 
  mutate(percent_alignment = percent_alignment / 100) %>% 
  mutate(million_reads_aligned = (M.Reads.After.Filtering..Trimmed.reads.)*percent_alignment) %>% 
  ggplot(., aes(x=Sample.ID, y=million_reads_aligned, fill=Treatment)) + geom_col() + facet_wrap(~Species) + theme_classic()
```
<img width="625" alt="Screen Shot 2023-08-14 at 4 51 09 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/48ce7fe6-33c4-47f9-a13b-9c2275db5487">

<img width="627" alt="Screen Shot 2023-08-14 at 4 51 23 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/0adc5f8f-1980-4f18-b1e7-8bc7315c0f5f">

Making data frame for reads aligned (to filter out anything less than 4 million reads)
```{r}
sequencing_data %>% 
  select(Species:Treatment,Raw.Reads..Million.,M.Reads.After.Filtering..Trimmed.reads.,Mapping.Rate..Post.trimming.alignment.) %>% 
  mutate(percent_alignment = gsub("%","", Mapping.Rate..Post.trimming.alignment.)) %>% 
  mutate(percent_alignment = as.numeric(percent_alignment)) %>% 
  mutate(percent_alignment = percent_alignment / 100) %>% 
  mutate(million_reads_aligned = (M.Reads.After.Filtering..Trimmed.reads.)*percent_alignment) -> alignment_rates

alignment_rates %>% 
  filter(million_reads_aligned > 4) %>%  # goes from 96 to 74 samples
  group_by(Species, Treatment, Genotype, Experiment.phase) %>% 
  summarise(n=n())
#this table shows the sample sizes of each Species + Treatment + Genotype combination at each time point

alignment_rates %>% 
  filter(million_reads_aligned > 4) %>%  # goes from 96 to 74 samples
  mutate(Experiment.phase = factor(Experiment.phase, levels = c("Pre-treatment", "last day of treatment"))) %>% 
  mutate(Genotype = factor(Genotype, c("A", "B", "C", "BC-8b", "MB-B", "SI-C"))) %>% 
  group_by(Species, Treatment, Genotype, Experiment.phase) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  ggplot(., aes(x=Genotype, y=n, fill=Treatment)) + 
  geom_col(position = position_dodge()) +
  theme_classic() +
  facet_grid(Species ~ Experiment.phase, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(labels=c("Control", "Variable"), values = c( "#60DBDB", "#F54A34"))
```

<img width="624" alt="Screen Shot 2023-08-14 at 4 51 56 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/46571088-4fdd-4b83-8645-c4e4a434e7e9">

## 2. DESeq2.Rmd

Then I made a separate R markdown file for the DESeq code.

```{r}
library(tidyverse)
library(DESeq2)
library(genefilter) #for filterfun and genefilter
library(factoextra) #for PCA and eigenvectors
library(vegan)
library(pheatmap)
```

Import and tidy gene counts matrix
```{r}
genecounts <- read.csv("../../results/Acer/gene_count_acerv_matrix.csv")
genecounts %>% column_to_rownames(var = "gene_id") -> genecounts

#tidy column names
colnames(genecounts) -> genecounts_colnames

gsub("^(.*?)_.*$", "\\1", genecounts_colnames) -> new_colnames

for (i in seq_along(new_colnames)) {
  colnames(genecounts)[colnames(genecounts) == colnames(genecounts)[i]] <- new_colnames[i]
}
```

import aligned reads sample metadata
```{r}
sequencing_data <- read.csv("../../RNA_extraction_sequencing_data.csv")

#obtaining list of Acer samples which passed read count
sequencing_data %>% 
  select(Species:Treatment,Sequence.File.Sample.Name,Raw.Reads..Million.,M.Reads.After.Filtering..Trimmed.reads.,Mapping.Rate..Post.trimming.alignment.) %>% 
  mutate(percent_alignment = gsub("%","", Mapping.Rate..Post.trimming.alignment.)) %>% 
  mutate(percent_alignment = as.numeric(percent_alignment)) %>% 
  mutate(percent_alignment = percent_alignment / 100) %>% 
  mutate(million_reads_aligned = (M.Reads.After.Filtering..Trimmed.reads.)*percent_alignment) %>% 
  filter(million_reads_aligned > 4) %>% 
  select(Species:Treatment,Sequence.File.Sample.Name,million_reads_aligned) %>% 
  filter(Species == "Acer") %>% 
  mutate(Sequence.File.Sample.Name = gsub("^(.*?)_.*$", "\\1", Sequence.File.Sample.Name)) %>% 
  mutate(Sequence.File.Sample.Name = gsub("-", ".", Sequence.File.Sample.Name)) -> Acer_samples_4M

#which ones need to be filtered out?
sequencing_data %>% 
  select(Species:Treatment,Sequence.File.Sample.Name,Raw.Reads..Million.,M.Reads.After.Filtering..Trimmed.reads.,Mapping.Rate..Post.trimming.alignment.) %>% 
  mutate(percent_alignment = gsub("%","", Mapping.Rate..Post.trimming.alignment.)) %>% 
  mutate(percent_alignment = as.numeric(percent_alignment)) %>% 
  mutate(percent_alignment = percent_alignment / 100) %>% 
  mutate(million_reads_aligned = (M.Reads.After.Filtering..Trimmed.reads.)*percent_alignment) %>% 
  filter(Species == "Acer") %>% 
  filter(million_reads_aligned < 4)
#Acer-072, Acer-089, Acer-108
```

<img width="985" alt="Screen Shot 2023-08-14 at 4 53 17 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/d9c64e47-8e02-4902-a547-676060d25492">

Check that there are no genes with 0 counts across all samples
```{r}
dim(genecounts) #33715    48
nrow(Acer_samples_4M) #45
#3 samples were filtered out

#Acer-072, Acer-089, Acer-108 from above

#manually remove them because you couldn't figure out how to use code to do it

samples_to_remove <- c("Acer.072", "Acer.089", "Acer.108")
genecounts_filt <- genecounts %>% select(-one_of(samples_to_remove))

#filter out counts that have zero
genecounts_filt <-genecounts_filt %>%
    mutate(Total = rowSums(.[, 1:45]))%>%
    filter(!Total==0)%>%
    dplyr::select(!Total)
nrow(genecounts_filt)

#went from 33715 to 14728 genes
```

Filter reads by proportion of samples containing cutoff value
```{r}
filt <- filterfun(pOverA(0.85,5)) #this means keep 85% of samples that have a count >5

tfil <- genefilter(genecounts_filt, filt)

keep <- genecounts_filt[tfil,]

gn.keep <- rownames(keep)

genecounts_filt <- as.data.frame(genecounts_filt[which(rownames(genecounts_filt) %in% gn.keep),])

#now at 4146 genes

#write.csv(genecounts_filt, "../../results/Acer/genecounts_filtered.csv")
```

Quality check of datasets to make sure row and column names match
```{r}
Acer_samples_4M %>% 
  column_to_rownames(var = "Sequence.File.Sample.Name") -> Acer_samples_4M

all(rownames(Acer_samples_4M) %in% colnames(genecounts_filt)) #TRUE 
all(rownames(Acer_samples_4M) == colnames(genecounts_filt)) # FALSE
```

Display order of metadata and gene count matrix, and reorder columns of genecounts_filt so they match the samples (need to do this before creating DESeq2 object because the rows and columns need to match EXACTLY)
```{r}
rownames(Acer_samples_4M)
colnames(genecounts_filt)

desired_order <- as.vector(rownames(Acer_samples_4M))

genecounts_filt %>% select(matches(desired_order)) -> genecounts_filt_ordered

all(rownames(Acer_samples_4M) == colnames(genecounts_filt_ordered)) #TRUE!
```

**Note:** I had to do a lot of data tidying to make sure the sample IDs matched across files.

set variables in metadata as factors
```{r}
Acer_samples_4M$Genotype <- factor(Acer_samples_4M$Genotype)
Acer_samples_4M$Treatment <- factor(Acer_samples_4M$Treatment)

#change genotype names to not have hyphens because DESeq2 doesn't like that
Acer_samples_4M$Genotype <- gsub("-", "_", Acer_samples_4M$Genotype)
#also need to get rid of spaces I think
Acer_samples_4M$Genotype <- gsub(" ", "", Acer_samples_4M$Genotype)

#make a new column for "Day_0" or "Day_29" to correspond to beginning and end sampling time points
Acer_samples_4M %>% 
  mutate(time_point = case_when(Experiment.phase == "Pre-treatment" ~ "Day_0",
                                Experiment.phase == "last day of treatment" ~ "Day_29")) -> Acer_samples_4M

Acer_samples_4M$time_point <- factor(Acer_samples_4M$time_point)
```

create matrix for DESeq
```{r}
data <- DESeqDataSetFromMatrix(countData = genecounts_filt_ordered, colData = Acer_samples_4M, design = ~ Genotype + Treatment*time_point)
```

This is the formula design I went with, because Treatment and time point are both things that I purposefully changed, while Genotype is more of a fixed effect (but something necessary to account for). 

I was worried about using time point as a factor instead of a continuous variable, but I looked at how [Dr. Kevin Wong did it in his 2019_Porites_rim_bleaching code](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/scripts/TagSeq/TagSeq_PCA_WGCNA.R) and he also included "Day" as a factor in his dds model.

estimate size factors
```{r}
SF.data <- estimateSizeFactors(data)
SF.data
print(sizeFactors(SF.data)) # everything is less than 4, so I can use vst
```

Apply variance stabilizing transformation to minimize effects of small counts and normalize wrt library size
```{r}
vst <- vst(data, blind = FALSE) #accounts for within group variability
head(assay(vst), 3)
```

scree plot for variance
```{r}
pca <- prcomp(t(assay(vst)))
fviz_eig(pca)
```

<img width="626" alt="Screen Shot 2023-08-14 at 4 56 16 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/0dda9328-f626-43c4-b01e-6b628eb6545b">

PCA
```{r}
plotPCA(vst, intgroup = c("Treatment"))
plotPCA(vst, intgroup = c("Genotype"))
plotPCA(vst, intgroup = c("time_point"))

vst_PCAdata <- plotPCA(vst, intgroup = c("Treatment", "Genotype", "time_point"), returnData = TRUE)
percentVar <- round(100*attr(vst_PCAdata, "percentVar")) 

#plot PCA of samples with all data
pca.centroids <- vst_PCAdata %>% 
  dplyr::select(Treatment, Genotype, time_point, PC1, PC2)%>%
  dplyr::group_by(Treatment, Genotype, time_point)%>%
  dplyr::summarise(PC1.mean = mean(PC1),
                   PC2.mean = mean(PC2))
find_hull <- function(vst_PCAdata) vst_PCAdata[chull(vst_PCAdata$PC1, vst_PCAdata$PC2), ]
hulls <- plyr::ddply(vst_PCAdata, "group", find_hull)

acerv_PCAplot <- ggplot(vst_PCAdata, aes(PC1, PC2, color=Genotype, shape=Treatment, alpha = time_point)) + 
   geom_point(size=3) +
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
   ggtitle("A. cervicornis - all genes") +
   theme_bw() + #Set background color
   theme(panel.border = element_blank(), # Set border
         axis.line = element_line(colour = "black"), #Set axes color
         plot.background=element_blank()) #Set the plot background
acerv_PCAplot
```

<img width="640" alt="Screen Shot 2023-08-14 at 4 56 39 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/46824fbd-abb8-43ae-b52d-8b1a8bc83e30">


This PCA is where things start getting hairy. There is an obvious separation of genotype, which is expected. There is also a huge time point separation, which is also expected but nice to see in the PCA so cleanly (although you can't really "see" it since I used opacity as the visualizing distinguisher since I ran out of other visuals). Treatment, the variable we most care about, is messy in there.

I wonder if it is appropriate/allowed to separate PCAs based on time point or genotype, and then see how treatment spreads out and whether it clusters? 

The only thing would be that if I were to continue down that route, I think I would still need to account for all the variables in the stats. But visualization-wise, I can separate things right? 

In any case I'm confused per usual. 
