---
layout: post
title: DESeq2 for Ch4 CCC samples (both location and genotype in dds model)
date: '2023-07-28'
categories: coding
tags: [coding, CCC_ch4]
---

I imported the quantified reads by using the multiqc general stats text file from the STAR alignment of trimmed samples and the readsPerGene.out.tab files generated from STAR alignment.

I followed [Natalia's code](https://github.com/China2302/SCTLD_RRC/blob/main/06a_gene_expression_analysis.Rmd) for this pipeline.

Here is the full R markdown I did to prep the count data (filtered for samples that had at least 4 million reads mapped):

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

reads_mapped <- read.delim("../multiqc_reports/multiqc_general_stats_STARalign_trimmed_updatedgfftake3.txt", sep = "\t")
head(reads_mapped)

reads_mapped<- reads_mapped %>% 
           dplyr::rename(percent_uniq = STAR_mqc.generalstats.star.uniquely_mapped_percent) %>%
  
               dplyr::rename(million_uniq = STAR_mqc.generalstats.star.uniquely_mapped) %>% 
  
               dplyr::rename(Sample_ID = Sample)

samples_selected_4M<- reads_mapped %>% filter(million_uniq > 4000000)
list_sampleIDs <- samples_selected_4M %>% select(Sample_ID) %>% as.list()

samples_selected_5M<- reads_mapped %>% filter(million_uniq > 5000000)
```

```{r}
# read in metadata
sample_metadata <- read.csv("../input_files/sample_metadata.csv")
names(sample_metadata)

#join read counts
sample_metadata_reads<- sample_metadata %>% left_join(samples_selected_4M, by="Sample_ID") 

sample_metadata_reads <- drop_na(sample_metadata_reads)
```

```{r}
# create a list of all files from samples that have at least 4M reads mapped

file_4M<- list.files("../star_trimmed_readsPerGeneoutTab_updatedgff3",
                  "*.gzReadsPerGene.out.tab", full.names = T)

countData_4M = data.frame(fread(file_4M[1]))[c(1,3)]

for(i in 2:length(file_4M)) {
        countData_4M = cbind(countData_4M, data.frame(fread(file_4M[i]))[3])
}

# Skip first 4 lines, count data starts on the 5th line
countData_4M = countData_4M[c(5:nrow(countData_4M)),]
colnames(countData_4M) = c("GeneID", gsub(paste0("../star_trimmed_readsPerGeneoutTab_updatedgff3"), "", file_4M))
colnames(countData_4M) = gsub("_trimmed_trimmed.fastq.gzReadsPerGene.out.tab", "",colnames(countData_4M))
colnames(countData_4M) = gsub("/", "ID_",colnames(countData_4M))
rownames(countData_4M) = countData_4M$GeneID
str(countData_4M) #33715 x 18

countData_4M = countData_4M[,c(2:ncol(countData_4M))] #this gets rid of the GeneID column but keeps GeneIDs as rownames 
str(countData_4M) #33715 x 17

write_rds(countData_4M, "countData_4M.rds")
```

Plot read counts
```{r}
ggplot(samples_selected_4M, aes(million_uniq)) + geom_histogram()

ggplot(samples_selected_4M, aes(percent_uniq)) + geom_histogram()

sample_metadata_reads %>% filter(!Genotype=="stag hybrid") %>% 
  ggplot(.) + 
  geom_bar(aes(x=Sample_ID, y=million_uniq, fill=Genotype), stat = "identity", position = position_dodge())

sample_metadata_reads %>% filter(!Genotype=="stag hybrid") %>% 
  ggplot(.) + 
  geom_bar(aes(x=Sample_ID, y=million_uniq, fill=Genotype), stat = "identity", position = position_dodge()) + 
  facet_wrap(~Location, scales = "free_x") +
  theme_classic() +
    theme(axis.text.x=element_text(angle=45,hjust=1))
```

<img width="459" alt="Screen Shot 2023-07-28 at 12 01 32 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/bc876e1e-1642-482e-b866-4712ae3a668a">

This figure is showing the number of uniquely mapped reads per sample, and is color coded by Genotype and separated by location. Unfortunately based on the poor quality of sequencing for some of the biological replicates, I don't think I can make a genotypic comparison for this (some genotypes only have one replicate at CCC). 

I wonder if I can remove genotype from the analysis and just focus on location?

Well for the first part of the DESeq2 code, I already made the DESeq2 object using both location and genotype. So here is the code for all of that (I followed [Kevin's](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/scripts/TagSeq/TagSeq_PCA_WGCNA.Rmd), [Ariana's](https://github.com/AHuffmyer/EarlyLifeHistory_Energetics/blob/master/Mcap2020/Scripts/TagSeq/Genome_V3/1_WGCNA_Mcap_V3.Rmd), [Jill's](https://github.com/JillAshey/SedimentStress/blob/master/RAnalysis/acerv/acerv_sub_DESeq2.Rmd), and [Natalia's](https://github.com/China2302/SCTLD_RRC/blob/main/06a_gene_expression_analysis.Rmd)):

Load libraries
```{r}
library(tidyverse)
library(DESeq2)
library(genefilter) #for filterfun and genefilter
library(factoextra) #for PCA and eigenvectors
library(vegan)
library(pheatmap)
```

Import data frames
```{r}
countData_4M<- readRDS("countData_4M.rds")
#remove ID 1096 because it is does not have any replicates
countData_4M %>% select(!ID_1096) -> countData_4M
dim(countData_4M) #33715 x 16

sample_metadata <- read_csv("../input_files/sample_metadata.csv")
sample_metadata$Sample_ID <- paste("ID_", sample_metadata$Sample_ID, sep="")
sample_metadata$Sample_ID <- gsub("_trimmed", "", sample_metadata$Sample_ID)

colnames_countData_4M<- as.data.frame(colnames(countData_4M)) %>% dplyr::rename('Sample_ID'='colnames(countData_4M)') %>% mutate(order = 1:16) 

#removing samples from metadata that did not meet the 4M cutoff
sample_metadata_4M <- left_join(colnames_countData_4M, sample_metadata, by = "Sample_ID") 

# Make row names to be SampleID 
rownames(sample_metadata_4M) <- sample_metadata_4M$Sample_ID

#change - to _ because the hyphen messes with DESeq2 somehow
sample_metadata_4M$Genotype <- gsub("-", "_", sample_metadata_4M$Genotype)
#also need to get rid of spaces I think
sample_metadata_4M$Genotype <- gsub(" ", "", sample_metadata_4M$Genotype)
```

Check that there are no genes with 0 counts across all samples
```{r}
nrow(countData_4M)
countData_4M_filt <-countData_4M %>%
    mutate(Total = rowSums(.[, 1:16]))%>%
    filter(!Total==0)%>%
    dplyr::select(!Total)
nrow(countData_4M_filt)

#went from 33715 to 19819 genes
```

Filter reads by proportion of samples containing cutoff value
```{r}
filt <- filterfun(pOverA(0.85,5)) #this means keep 85% of samples that have a count >5

tfil <- genefilter(countData_4M_filt, filt)

keep <- countData_4M_filt[tfil,]

gn.keep <- rownames(keep)

genecounts_filt <- as.data.frame(countData_4M_filt[which(rownames(countData_4M_filt) %in% gn.keep),])

#write.csv(genecounts_filt, "../results/genecounts_filtered.csv")
```

Quality check of datasets to make sure row and column names match
```{r}
all(rownames(sample_metadata_4M$Sample_ID) %in% colnames(genecounts_filt))
all(rownames(sample_metadata_4M$Sample_ID) == colnames(genecounts_filt))
#both return as TRUE
```

Display order of metadata and gene count matrix
```{r}
sample_metadata_4M$Sample_ID
colnames(genecounts_filt)
```

set location and genotype as factors
```{r}
sample_metadata_4M$Location <- factor(sample_metadata_4M$Location)
sample_metadata_4M$Genotype <- factor(sample_metadata_4M$Genotype)
```

### 1. DESeq accounting for both location and genotype ###

create matrix for DESeq
```{r}
data <- DESeqDataSetFromMatrix(countData = genecounts_filt, colData = sample_metadata_4M, design = ~ Location + Genotype)
```

Expression visualization

Text from Jill: "First we are going to log-transform the data using a variance stabilizing transforamtion (vst). This is only for visualization purposes. Essentially, this is roughly similar to putting the data on the log2 scale. It will deal with the sampling variability of low counts by calculating within-group variability (if blind=FALSE). Importantly, it does not use the design to remove variation in the data, and so can be used to examine if there may be any variability do to technical factors such as extraction batch effects. To do this we first need to calculate the size factors of our samples. This is a rough estimate of how many reads each sample contains compared to the others. In order to use VST (the faster log2 transforming process) to log-transform our data, the size factors need to be less than 4. Otherwise, there could be artefacts in our results."

My samples are all less than 4, so I'm good.
```{r}
SF.data <- estimateSizeFactors(data)
SF.data
print(sizeFactors(SF.data)) 

#Plot column sums according to size factor
plot(sizeFactors(SF.data), colSums(counts(SF.data)))
abline(lm(colSums(counts(SF.data)) ~ sizeFactors(SF.data) + 0))
#this is showing differences in sequencing depth
```

<img width="455" alt="Screen Shot 2023-07-28 at 12 24 46 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/f2669051-d0af-4e04-8777-4bebbe3a870a">

Apply variance stabilizing transformation to minimize effects of small counts and normalize wrt library size
```{r}
vst <- vst(data, blind = FALSE, fitType = 'local') #accounts for within group variability
head(assay(vst), 3)
```

1. Heatmap of sample-to-sample distances

```{r}
gsampleDists <- dist(t(assay(vst))) #calculate distance matix
gsampleDistMatrix <- as.matrix(gsampleDists) #distance matrix
rownames(gsampleDistMatrix) <- colnames(vst) #assign row names
colnames(gsampleDistMatrix) <- NULL #assign col names

pheatmap(gsampleDistMatrix, #plot matrix
         clustering_distance_rows=gsampleDists, #cluster rows
         clustering_distance_cols=gsampleDists) #cluster columns
```

<img width="449" alt="Screen Shot 2023-07-28 at 12 49 54 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/26cf6957-3f40-475b-bdc1-5e8013a9162f">
Idk what this plot means

2. Scree plot
```{r}
pca <- prcomp(t(assay(vst)))
fviz_eig(pca)
```
<img width="452" alt="Screen Shot 2023-07-28 at 12 53 38 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/617fb077-666a-441b-996e-4e1fd84dafdf">

3. PCA
```{r}
plotPCA(vst, intgroup = c("Location"))
plotPCA(vst, intgroup = c("Genotype"))

vst_PCAdata <- plotPCA(vst, intgroup = c("Location", "Genotype"), returnData = TRUE)
percentVar <- round(100*attr(vst_PCAdata, "percentVar")) #plot PCA of samples with all data
acerv_PCAplot <- ggplot(vst_PCAdata, aes(PC1, PC2, color=Location, shape=Genotype)) + 
   geom_point(size=3) +
   geom_text(aes(label=name),hjust=0, vjust=0) +
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
   ggtitle("A. cervicornis - all genes") +
   theme_bw() + #Set background color
   theme(panel.border = element_blank(), # Set border
         #panel.grid.major = element_blank(), #Set major gridlines
         #panel.grid.minor = element_blank(), #Set minor gridlines
         axis.line = element_line(colour = "black"), #Set axes color
         plot.background=element_blank()) #Set the plot background
acerv_PCAplot
```

<img width="451" alt="Screen Shot 2023-07-28 at 12 25 52 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/94a47f6f-338c-4f93-8f73-457144062f58">

This is the PCA for all genes, but what about the PCA for just DGEs? Jill separates out the DEGs at the end of her DESeq2 code and creates PCAs for just those. I'll try to follow the code step-by-step but I can't forget to do that at the end.

I also followed Kevin's code to try to add polygons, and this is what I got:
```{r}
pca.centroids <- vst_PCAdata %>% 
  dplyr::select(Location, Genotype, PC1, PC2)%>%
  dplyr::group_by(Location, Genotype)%>%
  dplyr::summarise(PC1.mean = mean(PC1),
                   PC2.mean = mean(PC2))
find_hull <- function(vst_PCAdata) vst_PCAdata[chull(vst_PCAdata$PC1, vst_PCAdata$PC2), ]
hulls <- plyr::ddply(vst_PCAdata, "group", find_hull)

acerv_PCAplot <- ggplot(vst_PCAdata, aes(PC1, PC2, color=Location, shape=Genotype)) + 
   geom_point(size=3) +
   geom_text(aes(label=name),hjust=0, vjust=0) +
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
   ggtitle("A. cervicornis - all genes") +
   theme_bw() + #Set background color
   theme(panel.border = element_blank(), # Set border
         #panel.grid.major = element_blank(), #Set major gridlines
         #panel.grid.minor = element_blank(), #Set minor gridlines
         axis.line = element_line(colour = "black"), #Set axes color
         plot.background=element_blank()) #Set the plot background
acerv_PCAplot

#add centroids
acerv_PCAplot + geom_polygon(data=hulls, alpha = 0.2, aes(color = Location, fill = Location, lty = Genotype)) +
  geom_point(aes(x=PC1.mean, y=PC2.mean,color=Location, shape = Genotype), data=pca.centroids, size=4, show.legend=FALSE) + 
  scale_linetype_manual(values = c("solid", "dashed", "dotted"))
  # scale_shape_manual(values=c(15, 17, 19)) +
  # theme(legend.text = element_text(size=8), 
  #       legend.position=c(0.95,0.85),
  #       plot.background = element_blank(),
  #       legend.title = element_text(size=10), 
  #       plot.margin = margin(1, 1, 1, 1, "cm"),
  #       axis.text = element_text(size=18), 
  #       title = element_text(size=25, face="bold"), 
  #       axis.title = element_text(size=18))
```
<img width="452" alt="Screen Shot 2023-07-28 at 1 15 23 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/bf986686-e2b2-4392-a967-bfcab44b4a2d">

Which at first glance looks good, but then when you look to the far right you see a tiny triangle and I think it is trying to create polygons based on the Genotype:Location combinations (i.e. MiamiBeach_C at CCC is the large red triangle). I want to create polygons for just location, because I think that's the important separator, but when I try just removing genotype from the code, I get straight lines instead. So idk how to do that/if it's even meaningful because Kevin uses this code to then compare shifts of expression patterns over time. 

He then runs a PERMANOVA:
```{r}
test<-t(assay(vst))
test<-as.data.frame(test)

test$Sample_ID <-rownames(test)
test$Location <- sample_metadata_4M$Location[match(test$Sample_ID, sample_metadata_4M$Sample_ID)]
test$Genotype <- sample_metadata_4M$Genotype[match(test$Sample_ID, sample_metadata_4M$Sample_ID)]
```

```{r}
dim(test)
scaled_test <-prcomp(test[c(1:1304)], scale=TRUE, center=TRUE)
fviz_eig(scaled_test)
vegan <- scale(test[c(1:1304)])

permanova<-adonis2(vegan ~ Location*Genotype, data = test, method = 'eu')
```
Result: not significant
<img width="461" alt="Screen Shot 2023-07-28 at 1 24 58 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/dc27ea86-398a-4482-a545-f27a5c50ab9c">

But then the number of replicates per Location:Genotype is very low, only 1 in some cases. So I don't think this statistical test is even meaningful. After I run the DESeq I will redo this whole part with just Location as a variable.

Run DESeq
```{r}
DEG_locgen <- DESeq(data, fitType = 'local')
DEG_locgen_res <- results(DEG_locgen, alpha = 0.05)
summary(DEG_locgen_res)
```
Results:
<img width="357" alt="Screen Shot 2023-07-28 at 1 28 14 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/c9411f4a-30ab-43fe-a8ac-6e3f1d9129a1">

L. O. L. 

But when I specify the contrast of interest (Location nursery vs. CCC), I get:
```{r}
resultsNames(DEG_locgen)
[1] "Intercept"                         
[2] "Location_nursery_vs_CCC"           
[3] "Genotype_MiamiBeach_C_vs_Cheetos_B"
[4] "Genotype_SunnyIsles_E_vs_Cheetos_B"

DEG_Nursery_vs_CCC <- results(DEG_locgen, contrast = c("Location", "nursery", "CCC"), alpha = 0.05)
summary(DEG_Nursery_vs_CCC)
```
<img width="359" alt="Screen Shot 2023-07-28 at 1 34 11 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/cc10cc02-a311-4d4c-97c3-c4f1ee21adcf">

I get some numbers so that's good.

Compare Nursery vs. CCC
```{r}
DEG_Nursery_vs_CCC <- as.data.frame(DEG_Nursery_vs_CCC)
DEG_Nursery_vs_CCC["Location_Compare"] <- "NurseryvsCCC"
#write.csv(DEG_Nursery_vs_CCC, file = "Nursery_vs_CCC_allgenes_201230728.csv")

DEG_Nursery_vs_CCC.sig.num <- sum(DEG_Nursery_vs_CCC$padj<0.05, na.rm=T)
#455 DEGs

#get a list of just those genes
DEG_Nursery_vs_CCC.sig <- subset(DEG_Nursery_vs_CCC, padj <0.05) # identify and subset significant pvalues
DEG_Nursery_vs_CCC.sig["Location_Compare"] <- "NurseryvsCCC" # adding treatment comparison column
DEG_Nursery_vs_CCC.sig.list <- data[which(rownames(data) %in% rownames(DEG_Nursery_vs_CCC.sig)),] # subset list of significant genes from original count data 
DEG_Nursery_vs_CCC.sig.list <- as.data.frame(counts(DEG_Nursery_vs_CCC.sig.list)) # make list of sig gene counts into a df
DEG_Nursery_vs_CCC.sig.list_full <- cbind(DEG_Nursery_vs_CCC.sig, DEG_Nursery_vs_CCC.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_Nursery_vs_CCC.sig.list_full, file = "DEG_Nursery_vs_CCC.sig.list_full_20230728.csv") # write out csv
```

Variance stabilized transformation for just DEGs
```{r}
DEG_Nursery_vs_CCC.sig.list <- data[which(rownames(data) %in% rownames(DEG_Nursery_vs_CCC.sig.list)),] 
# turn back into formal class DESeqTransform or else vst will not run

SFtest <- estimateSizeFactors(DEG_Nursery_vs_CCC.sig.list)
print(sizeFactors(SFtest)) #everything is less than 4 so we can do vst

DEG_Nursery_vs_CCC.sig.vst <- varianceStabilizingTransformation(DEG_Nursery_vs_CCC.sig.list, blind = FALSE, fitType = 'local') # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
```

PCA plot of DEGs
```{r}
acerv_sub_DEG_PCA <- plotPCA(DEG_Nursery_vs_CCC.sig.vst, intgroup = c("Location"), returnData=TRUE)
percentVar_pca_acerv_sub <- round(100*attr(acerv_sub_DEG_PCA, "percentVar")) #plot PCA of samples with all data

acerv_sub_DEG_PCA_plot <- ggplot(acerv_sub_DEG_PCA, aes(PC1, PC2, color = Location)) +
  geom_point(size=6) +
  xlab(paste0("PC1: ",percentVar_pca_acerv_sub[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_acerv_sub[2],"% variance")) +
  #coord_fixed() +
  ggtitle(label = "A. cervicornis") +
  theme_bw() + #Set background color
  theme(legend.text = element_text(size=8), 
        #legend.position="none",
        plot.background = element_blank(),
        #legend.title = element_text(size=18, face="bold"), 
        legend.title=element_blank(),
        axis.text = element_text(size=8), 
        axis.title = element_text(size=10,  face="bold"), 
        axis.title.y = element_text(vjust=-1.5),
        plot.title = element_text(size = 15, face = "italic", hjust = 0.5))
acerv_sub_DEG_PCA_plot # PCA plot is of differentially expressed genes only
```
<img width="453" alt="Screen Shot 2023-07-28 at 2 46 56 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/8fe078da-c486-4c93-94cc-e5af2e0f8c28">

2. Dispersion Plot
```{r}
plotDispEsts(DEG_locgen, main="Dispersion plot")
```
<img width="450" alt="Screen Shot 2023-07-28 at 2 48 15 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/d7ace87e-9837-4d56-bd81-4b7045127e86">

3. Cook's Distance
```{r}
boxplot(log10(assays(DEG_locgen)[["cooks"]]), range=0, las=0, main="Cook's distance")
```
<img width="449" alt="Screen Shot 2023-07-28 at 2 49 05 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/abf97cb3-6e67-4b3e-9535-174f60927cbf">


