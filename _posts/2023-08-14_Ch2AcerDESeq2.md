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

**08.15.2023 update:** 

Let's try separating the two time points into their own PCAs so it is easier to distinguish genotype x treatment. 

filter for just the first time point
```{r}
Acer_samples_4M %>% filter(time_point == "Day_0") -> day_0_samples

list_day0 <- rownames(day_0_samples)

genecounts_filt_ordered %>% select(matches(list_day0)) -> day_0_genecounts
```

create matrix for DESeq
```{r}
day0_dds <- DESeqDataSetFromMatrix(countData = day_0_genecounts, colData = day_0_samples, design = ~ Treatment + Genotype)
```

estimate size factors
```{r}
SF_day0 <- estimateSizeFactors(day0_dds)
print(sizeFactors(SF_day0)) # everything is less than 4, so I can use vst
```


Apply variance stabilizing transformation to minimize effects of small counts and normalize wrt library size
```{r}
vst_day0 <- vst(day0_dds, blind = FALSE) #accounts for within group variability
```

scree plot for variance
```{r}
pca_day0 <- prcomp(t(assay(vst_day0)))
fviz_eig(pca_day0)
```

PCA
```{r}
plotPCA(vst_day0, intgroup = c("Treatment"))
plotPCA(vst_day0, intgroup = c("Genotype"))

vst_day0_PCAdata <- plotPCA(vst_day0, intgroup = c("Treatment", "Genotype"), returnData = TRUE)
percentVar <- round(100*attr(vst_day0_PCAdata, "percentVar")) 

#plot PCA of samples with all data
pca.centroids_day0 <- vst_day0_PCAdata %>% 
  dplyr::select(Treatment, Genotype, PC1, PC2)%>%
  dplyr::group_by(Treatment, Genotype)%>%
  dplyr::summarise(PC1.mean = mean(PC1),
                   PC2.mean = mean(PC2))
find_hull_day0 <- function(vst_day0_PCAdata) vst_day0_PCAdata[chull(vst_day0_PCAdata$PC1, vst_day0_PCAdata$PC2), ]
hulls_day0 <- plyr::ddply(vst_day0_PCAdata, "group", find_hull_day0)

ggplot(vst_day0_PCAdata, aes(PC1, PC2, color=Genotype, shape=Treatment)) + 
   geom_point(size=3) +
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
   ggtitle("Day 0 Acer - all genes") +
   theme_bw() + #Set background color
   theme(panel.border = element_blank(), # Set border
         axis.line = element_line(colour = "black"), #Set axes color
         plot.background=element_blank()) #Set the plot background
```

<img width="504" alt="Screen Shot 2023-08-15 at 10 24 50 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/8f772705-97a5-4a1d-a559-7b2ccfb1b605">

Let's try to make polygons around the treatments specifically.

So if I try to make the polygons around the treatments, it looks weird because there is a huge separation of genotype and it's trying to connect all those dots:

<img width="814" alt="Screen Shot 2023-08-15 at 10 35 59 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/8179d36f-6fd3-4f24-91db-fcd2dc0ed3d7">

We don't want that. I'll try to flip it so it groups by genotype instead:

<img width="804" alt="Screen Shot 2023-08-15 at 10 38 26 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/1b84baa1-7096-4da1-a611-89ebaa29916b">

I don't get what this is doing either. It also looks weird?

There is so much clustering of genotype. 

Ok this is kind of getting to where I was envisioning, where the polygons are around the subgroupings. It just looks chaotic and messy.

<img width="925" alt="Screen Shot 2023-08-15 at 10 48 42 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/5b1b100b-af5b-4a56-831f-0fcbe729dd98">

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

ggplot(vst_PCAdata, aes(PC1, PC2, color=Treatment, shape=Treatment)) + 
   geom_point(size=3) +
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
   ggtitle("A. cervicornis - all genes") +
   theme_bw() + #Set background color
   theme(panel.border = element_blank(), # Set border
         axis.line = element_line(colour = "black"), #Set axes color
         plot.background=element_blank()) + #Set the plot background
  geom_polygon(data = hulls, alpha = 0.2, aes(color = group, fill = group))
```

I might need to look more into what the hulls code is actually calculating, and that could make the polygons cleaner. 

What about that ellipse code?

<img width="566" alt="Screen Shot 2023-08-15 at 10 56 24 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/ae598579-109d-4b0f-afd2-9e48902f9381">

```{r}
ggplot(vst_PCAdata, aes(PC1, PC2, color=Genotype, shape=Treatment)) + 
   geom_point(size=3) +
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
   ggtitle("A. cervicornis - all genes") +
   theme_bw() + #Set background color
   theme(panel.border = element_blank(), # Set border
         axis.line = element_line(colour = "black"), #Set axes color
         plot.background=element_blank()) + 
  stat_ellipse(aes(PC1, PC2, group=Genotype), type = "norm")
```

<img width="567" alt="Screen Shot 2023-08-15 at 10 58 08 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/fc719563-4a64-4622-bb7b-60173ac40968">


```{r}
ggplot(vst_PCAdata, aes(PC1, PC2, color=Treatment, shape=Genotype)) + 
   geom_point(size=3) +
  scale_color_manual(labels=c("Control", "Variable"), values = c( "#60DBDB", "#F54A34")) +
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
   ggtitle("A. cervicornis - all genes") +
   theme_bw() + #Set background color
   theme(panel.border = element_blank(), # Set border
         axis.line = element_line(colour = "black"), #Set axes color
         plot.background=element_blank()) + 
  stat_ellipse(aes(PC1, PC2, group=time_point, lty = time_point), type = "norm")
```

<img width="567" alt="Screen Shot 2023-08-15 at 10 59 34 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/b10c7c2f-d7d8-466e-b82d-3d22b953f5f0">

```{r}
ggplot(vst_PCAdata, aes(PC1, PC2, color=Genotype, shape=time_point)) + 
   geom_point(size=3) +
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
   ggtitle("A. cervicornis - all genes") +
   theme_bw() + #Set background color
   theme(panel.border = element_blank(), # Set border
         axis.line = element_line(colour = "black"), #Set axes color
         plot.background=element_blank()) + 
  stat_ellipse(aes(PC1, PC2, group=Treatment, lty = Treatment), type = "norm")
```

## Building PERMANOVA model

Conduct PERMANOVA
```{r}
test<-t(assay(vst)) #remember that formula is ~ Treatment*time_point + Genotype
test<-as.data.frame(test)

test$Sample_ID <-rownames(test)
test$Treatment <- Acer_samples_4M$Treatment[match(test$Sample_ID, rownames(Acer_samples_4M))]
test$Genotype <- Acer_samples_4M$Genotype[match(test$Sample_ID, rownames(Acer_samples_4M))]
test$time_point <- Acer_samples_4M$time_point[match(test$Sample_ID, rownames(Acer_samples_4M))]
```

Build PERMANOVA model
```{r}
dim(test) #45 4150
scaled_test <-prcomp(test[c(1:4146)], scale=TRUE, center=TRUE) #subtract 4 from the dim because you made 4 columns for metadata
fviz_eig(scaled_test)
vegan <- scale(test[c(1:4146)])

permanova<-adonis2(vegan ~ Treatment*time_point + Genotype, data = test, method = 'eu')

print(permanova) #everything is significant woooooo
```

<img width="453" alt="Screen Shot 2023-08-15 at 12 17 11 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/a0c651b4-e31a-453b-ace7-2c27035845a8">


## Running DESeq2

un DESeq2
```{r}
DEG_all <- DESeq(data)
DEG_all_res <- results(DEG_all, alpha = 0.05)
summary(DEG_all_res)

resultsNames(DEG_all)

DEG_treatment <- results(DEG_all, contrast = c("Treatment", "variable", "control"), alpha = 0.05)
summary(DEG_treatment)

DEG_timepoint <- results(DEG_all, contrast = c("time_point", "Day_29", "Day_0"), alpha = 0.05)
summary(DEG_timepoint)
```
DEG_treatment
<img width="355" alt="Screen Shot 2023-08-15 at 12 18 09 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/f54683c1-d965-4596-9d3d-6ba094bb32c6">

DEG_timepoint
<img width="393" alt="Screen Shot 2023-08-15 at 12 18 27 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/26ce9c33-09a0-41df-84e5-c06b0ed346b0">


I don't think these results are informative. I want a result that says "variable versus control at day 29 (so after the treatment to get an effect of treatment?)". or do I want to also see what they looked like at day 0? or do I do an LRT where I use day 0 as the "baseline" and then compare the effect of treatment after the 29 days. I don't necessarily care about the differential gene expression on day 0. That would be the sort of "baseline" gene expression differences between genotypes.

Let's try running LRT to see if the reduced model of treatment + time_point is versus the full model is significant? (essentially accounting for the baseline differences of both at the start of the experiment that are not due to treatment)

DESeq2 with LRT
```{r}
data <- DESeqDataSetFromMatrix(countData = genecounts_filt_ordered, colData = Acer_samples_4M, design = ~ Treatment*time_point + Genotype)

ddsLRT <- DESeq(data, test="LRT", reduced = ~Treatment + time_point)

 #this reduced function is including treatment + time_point (but not the interaction), which means than any differences at time point 0 will be controlled for. (see Michael Love's response on this thread: https://support.bioconductor.org/p/62684/  "Using a reduced formula of ~ time + treat means that genes which show a consistent difference from time 0 onward will not have a small p value. This is what we think you want, because differences between groups at time 0 should be controlled for, e.g. random differences between the individuals chosen for each treatment group which are observable before the treatment takes effect.")

resLRT <- results(ddsLRT)
summary(resLRT, alpha = 0.05)

resultsNames(ddsLRT)

#identifying significant genes
resSig_LRT<-resLRT[which(resLRT$padj<0.05), ]
```

<img width="368" alt="Screen Shot 2023-08-15 at 12 19 57 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/f427ae0a-7b43-4bf6-9f88-ea535ae8f6dd">


Should I proceed with this model? I don't think the WGCNA part needs the LRT vs. Wald test distinction. And when I do the vst transformation, that is before I run the DESeq2 anyways so that shouldn't affect the PERMANOVA either. Or should I be running the PERMANOVA only on the significant DGEs????????? 

**3:43 pm update**:

So I looked through Kevin Wong's [Porites_Rim_Bleaching](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/scripts/TagSeq/TagSeq_PCA_WGCNA.Rmd) code and he didn't do a PERMANOVA on his DESeq2 stuff, so I'm not sure where I got that from to be honest. I think he did a PERMANOVA on physiological traits... [see his Rmd for the physiology analysis and PERMANOVA starting at line 1113](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/scripts/Physiology/Physiology_Analysis.Rmd). 

Ok so I see what I did here. He had made a PCA code to make centroids and generally display all his physiology data as a PCA. Then he ran a PERMANOVA on that. I don't think I can do a PERMANOVA on differential gene expression analysis... even though I'm making a PCA? Or wait did I talk to him about this a long time ago and he said it was fine because a PERMANOVA is non-parametric?

Yes, I think that's what it was, because when I presented the wound healing data as a PCA he asked what the stats of that was. So I then [added PERMANOVAs for them here but they were not significant](https://github.com/ademerlis/woundhealingPdam/blob/main/R_markdowns/PCA_permanovas_individ_hours.Rmd). 

Ok circling back to the point of this. The goal is to create a *meaningful* differential expression analysis. I don't think the current formula I'm using is perfectly encapturing all the interactions I want to look at -- my reason for believing this is when I run `resultsNames(dds)` for the original dds object and the dds_LRT, I get the same names (which don't look right):

<img width="279" alt="Screen Shot 2023-08-15 at 3 48 52 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/9bd604ec-e4b2-47f5-8e33-411427e76609">

What I want it to look like is something like:

control_day29 versus variable_day29
genotype1_variable_day0 versus genotype1_variable_day29
variable_day0 versus variable_day29 

I would want the differences between control_day0 and variable_day0 to be accounted for -- which is I think why I use the LRT (and have that as the "reduced" factor). 

I was looking for some DESeq code from [Dr. Carly Kenkel's lab group](https://dornsife.usc.edu/labs/carlslab/data/) because I know she does coral transcriptomics work, and I came across PhD student [Yingqi Zhang's GitHub repo](https://github.com/yingqizhang/OfavLarvae/blob/main/DESeq2/DESeq2.R) for transcriptomics of Ofav larvae. 

In that, she does something that I hadn't thought about doing before for a DESeq2 model, but that makes complete sense:

```{r}
# combine the factors of interest into a single factor with all combinations of the original factors
dds4 <- DESeqDataSetFromMatrix(countData = cts_ms[,1:46],colData = colData_ms,design = ~ 1)
dds4$group <- factor(paste0(dds4$Trmt, dds4$Type))

# change the design to include just this factor, e.g. ~ group
design(dds4) <- ~ group
dds4 <- DESeq(dds4)
res4=results(dds4)
summary(res4) # 405 up, 301 down, 7571 low count

res_inshoretrmt <- results(dds4, contrast=c("group","HeatInshore","CtrlInshore"))
summary(res_inshoretrmt) # 133 up, 144 down
```

So what I need to do is make a matrix that has like day_treatment_genotype and every combination of that. It looks like based on her code that you add that to the design after the fact.

```{r}
# combine the factors of interest into a single factor with all combinations of the original factors
dds_combo <- DESeqDataSetFromMatrix(countData = genecounts_filt_ordered,colData = Acer_samples_4M,design = ~ 1)
dds_combo$group <- factor(paste0(dds_combo$Treatment, "_", dds_combo$time_point, "_", dds_combo$Genotype)) #12 levels

# change the design to include just this factor, e.g. ~ group
design(dds_combo) <- ~ group

dds_combo <- DESeq(dds_combo) #this is a wald test though of doing pairwise comparisons for all combos

results_combo <- results(dds_combo)
summary(results_combo) #427 up, 384 down
```

<img width="349" alt="Screen Shot 2023-08-15 at 4 01 42 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/dcfbd628-54be-4e20-84f7-d7c18c97698b">


So that works and is great, but what about the "reduced" formula with LRT? I'm not sure how to specify that when genotype is in the formula. because I want it to be reduced by treatment + day_0 essentially. maybe I filter out everything with day 0? but then is that even informative? so maybe I don't need to do an LRT for this... unless I had an interaction term that was treatment x time point. 

<img width="464" alt="Screen Shot 2023-08-15 at 4 07 54 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/f6f140f2-f4dd-4859-a224-798c729b862f">

With these comparisons, I can't see doing a reduced model being informative. 

Also, wait these are not all the possible comparisons that can be made. It is actually every other comparison versus the first one, which is "control_day0_BC-8b". So I'm confused actually how this design was made. 

Maybe this doesn't work.

Ok wait I just played around with the formula and combinations, and this might have worked.

```{r}
# Trying to create a combination where genotype is not part of each comparison, just the main factors of treatment and time point are created as a combo factor so I can do every comparison of those combinations

Acer_samples_4M$Genotype <- as.factor(Acer_samples_4M$Genotype)

dds_combo <- DESeqDataSetFromMatrix(countData = genecounts_filt_ordered,colData = Acer_samples_4M,design = ~ 1)
dds_combo$group <- factor(paste0(dds_combo$Treatment, "_", dds_combo$time_point)) #12 levels

# change the design to include just this factor, e.g. ~ group
design(dds_combo) <- ~ (group + Genotype + group:Genotype)

dds_combo <- DESeq(dds_combo) #this is a wald test though of doing pairwise comparisons for all combos

results_combo <- results(dds_combo)
summary(results_combo) #9 up, 0 down

resultsNames(dds_combo)
```

<img width="585" alt="Screen Shot 2023-08-15 at 4 40 08 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/56d539ce-dcc7-420b-8a95-2aa635b0e64d">

So that looks kind of like what I was envisioning, except I'm not really sure what the genotypes are doing. (like what does "groupcontrol_Day_29.GenotypeMB_B" mean? the interaction of untreated_day29 and genotype MB-B?). 

But something preliminary that is exciting (but seems too good to be true). When I run the comparison of control day 0 versus variable day 0, i get NO DGES!!!! which is insane? maybe if i were to do within genotype comparisons that wouldn't be the case and would make more sense. 

<img width="620" alt="Screen Shot 2023-08-15 at 4 42 17 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/7902423c-347d-4b4c-a493-2557d349fb4d">


But when I do variable vs control on day 29, I get a bunch of DGEs!

<img width="593" alt="Screen Shot 2023-08-15 at 4 43 14 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/348a09fd-0ce3-489c-8227-719deffa79d6">

So that's something anyways.
