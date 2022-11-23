---
layout: post
title: DEGPatterns Gene Clustering 
date: '2022-11-22'
categories: Coding
tags: [Coding]
---

I have been trying to recreate this figure but I'm having a lot of trouble getting it.

![Cluster_Plot](https://user-images.githubusercontent.com/56000927/203420486-0675395f-058a-4221-aae7-0a450c92c450.png)

### The code used is below:

```{r}

library(knitr)
library(DESeq2) #needs BiocManager
library(DEGreport) #BiocManager
library(edgeR) #BiocManager
library(limma) #BiocManager
library(Biobase) 
library(tidyverse) 
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(apeglm) #BiocManager
library(topGO) #BiocManager
library(GO.db) #BiocManager
library(Rgraphviz) #BiocManager
library(genefilter)
library(ggpubr)

#citation("tximport")

#reading in txt file as a table with the headers as the header in the original file
countsmatrix <-read.table("WH_Pdam.Rmatrix.txt",header=TRUE)

#removing the Chromosome column from the dataset, and changing the rows to be the gene names
countsmatrix <-countsmatrix %>% 
  dplyr::select(-c("Chr")) %>% 
  column_to_rownames("Geneid")

#removing hour 5 from the dataset because it only has 1 wounded coral, so it is not comparable 
countsmatrix <- countsmatrix %>% 
  dplyr::select(!X5331C1:X5335Y)

#function found from stack overflow
#tidying column names so they don't include an X anymore (R puts an X in front of columns that start with numbers)
destroyX = function(es) {
  f = es
  for (col in c(1:ncol(f))){ #for each column in dataframe
    if (startsWith(colnames(f)[col], "X") == TRUE)  { #if starts with 'X' ..
      colnames(f)[col] <- substr(colnames(f)[col], 2, 100) #get rid of it
    }
  }
  assign(deparse(substitute(es)), f, inherits = TRUE) #assign corrected data to original name
}

destroyX(countsmatrix)

#making the metadata file for the sample and treatment information
metadata = data.frame(sample=colnames(countsmatrix),
                condition = stringr::str_detect(pattern = ".*C.*",string = colnames(countsmatrix)),
                hour = stringr::str_replace(pattern = "(.).*",replacement="\\1",string = colnames(countsmatrix)),
                id = stringr::str_replace(pattern="(...).*",replacement="\\1",string=colnames(countsmatrix)))

#changing TRUE and FALSE for condition to Control and Wounded
metadata$condition[str_detect(metadata$condition,"TRUE")] <- "Control"
metadata$condition[str_detect(metadata$condition,"FALSE")] <- "Wounded"

metadata <- metadata %>% column_to_rownames("sample")

genefeatures <- read.delim(file = "pdam_genome_IDInfo.gff", header = F)
head(genefeatures)
colnames(genefeatures) <- c("IDGeneInfo")
rownames(genefeatures) <- rownames(countsmatrix)

#removing gene name from the second column
gene_name_split <- str_split_fixed(genefeatures$IDGeneInfo, " ", 2)
print(head(gene_name_split))
colnames(gene_name_split)<- c("Gene_ID","Gene_Function")
head(gene_name_split)
gene_name_split_df<-data.frame(gene_name_split)
head(gene_name_split_df)
#convert the Gene_ID column to the row headers
gene_name_split_df %>%
     remove_rownames() %>%
 column_to_rownames(var = 'Gene_ID')-> Pdam_Gene_Names_Info 
#%>%
  #write.csv(file="Pdam_Gene_Names_Info.csv")

#before creating your dds variable, you must ensure that the rownames of your metadata file are the same (and in the same order) as the columns for your counts matrix 
all(rownames(metadata) == colnames(countsmatrix))

#DESeq2 for Wald test from a counts matrix; requires count data, metadata, and the experimental design
dds=DESeq2::DESeqDataSetFromMatrix(countData = countsmatrix, colData = metadata, design = ~condition+hour+condition:hour)
dds <- estimateSizeFactors(dds)

#prefiltering recommended by DESeq2 Vignette (removes anything with reads below 10)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Run DESeq2
dds = 
dds <- readRDS(file = "WH_LRT_DESeq.rds")
res <- results(dds)
res
summary(res)


#DESeq2 variable for Wald test
dds_Wald<-DESeq2::DESeq(dds)

#DESeq2 variable for LRT test
dds_LRT<-DESeq2::DESeq(dds,reduced=~hour,test="LRT")

#LRT is the Likelihood Ratio Test, used to identify differentially expressed genes (DEGs) between control vs. treatment over the 5 time points (hours 0-5)

#LRTs are very "generous", so a more stringent alpha value than 0.05 is needed
res01<- results(dds_LRT, alpha = 0.01)

#annotate with Gene Function
results_LRT_table<-merge(as.data.frame(res01),Pdam_Gene_Names_Info,by='row.names',all=FALSE)%>%
  column_to_rownames("Row.names")

all(rownames(results_LRT_table) == rownames(Pdam_Gene_Names_Info))
#FALSE - not all genes are in res01

write.csv(results_LRT_table,file="Results_LRT_Table.csv")

#filtering by a significance value of p < 0.01 for gene expression 
results_table_LRT_01 <- results_LRT_table %>% rownames_to_column("Geneid")%>% filter(padj < 0.01)

# Get sig gene lists
sigLRT_genes <- results_table_LRT_01 %>% pull(Geneid)
length(sigLRT_genes)
#67 

write.csv(results_table_LRT_01,file="Results_Table_LRT_01.csv")

#save LRT DEseq as an object
saveRDS(dds_LRT, "WH_LRT_DESeq.rds")

### Compute pairwise correlation values
dds_vst<- vst(dds_LRT,blind=FALSE)

rld<-rlog(dds_LRT,blind=TRUE)
rld_mat <- assay(rld)
# Input is a matrix of log transformed values

rld_cor <- cor(rld_mat)    ## cor() is a base R function

head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames

#DestoryX formula for row
destroyX_row = function(es) {
  f = es
  for (row in c(1:nrow(f))){ #for each column in dataframe
    if (startsWith(rownames(f)[row], "X") == TRUE)  { #if starts with 'X' ..
      rownames(f)[row] <- substr(rownames(f)[row], 2, 100) #get rid of it
    }
  }
  assign(deparse(substitute(es)), f, inherits = TRUE) #assign corrected data to original name
}

rld_cor<-destroyX_row(rld_cor)

### Plot heatmap - need to figure out how to annotate for condition
heat.colors <- brewer.pal(6, "Blues")
pheatmap(rld_cor, color = heat.colors, border_color=NA, fontsize = 10, 
  		fontsize_row = 10, height=10, width=10, filename = "Heatmap_VST_DDS.jpeg")

#need to change 'hour' in metadata from chr to number 
metadata$hour<-as.factor(metadata$hour)

clustering_sig_genes <- results_table_LRT_01 %>%
                  arrange(padj) %>%
                  head(n=1000)

cluster_rlog <- rld_mat[clustering_sig_genes$Geneid, ]
destroyX(cluster_rlog)

# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups 

clusters <- degPatterns(cluster_rlog, metadata, time="hour", col="condition") #warnings can be ignored

#editing the plot for publication
cluster_plot<-clusters$plot+
  xlab(paste0("Hour"))+
  theme(text=element_text(size=12))+
  theme(legend.key.size=unit(0.5,"cm"))+
  theme(legend.title=element_text(size=12))+
  labs(fill="Condition")+
  scale_color_discrete(name="Condition")+
  theme_classic(base_size=12,base_family="serif")
ggsave(filename="Cluster_Plot.png")
  
# What type of data structure is the `clusters` output?
class(clusters)

# Let's see what is stored in the `df` component
head(clusters$df)

#extract raw data from clusters
clustersdata<-as.data.frame(clusters$raw)

#cluster 1 table
cluster1<-clustersdata[clustersdata$cluster=='1',]
#clean up table and annotate with gene function
cluster1<-subset(cluster1,merge=="Control0",select=-c(merge,condition,cluster,id,hour))
#picked Control0 arbitrarily (cluster doesn't differ by condition or hour)
rownames(cluster1) <- cluster1$genes 
cluster1<-merge(cluster1,Pdam_Gene_Names_Info,by='row.names',all=TRUE)%>%
  na.omit(cluster1)
cluster1<-subset(cluster1,select=-c(Row.names,value))
#merge with cluster_sig_genes table to sort by significance
colnames(cluster1)[colnames(cluster1)=="genes"]<-"Geneid"
cluster1<-merge(cluster1,clustering_sig_genes,by='Geneid')
cluster1<-subset(cluster1,select=-c(Gene_Function.x))

#cluster 2 table
cluster2<-clustersdata[clustersdata$cluster=='2',]
#clean up table and annotate with gene function
cluster2<-subset(cluster2,merge=="Control0",select=-c(merge,condition,cluster,id,hour))
rownames(cluster2) <- cluster2$genes 
cluster2<-merge(cluster2,Pdam_Gene_Names_Info,by='row.names',all=TRUE)%>%
  na.omit(cluster2)
cluster2<-subset(cluster2,select=-c(Row.names,value))
#merge with cluster_sig_genes
colnames(cluster2)[colnames(cluster2)=="genes"]<-"Geneid"
cluster2<-merge(cluster2,clustering_sig_genes,by='Geneid')
cluster2<-subset(cluster2,select=-c(Gene_Function.x))

#cluster 3 table
cluster3<-clustersdata[clustersdata$cluster=='3',]
#clean up table and annotate with gene function
cluster3<-subset(cluster3,merge=="Control0",select=-c(merge,condition,cluster,id,hour))
rownames(cluster3) <- cluster3$genes 
cluster3<-merge(cluster3,Pdam_Gene_Names_Info,by='row.names',all=TRUE)%>%
  na.omit(cluster3)
cluster3<-subset(cluster3,select=-c(Row.names,value))
#merge with cluster_sig_genes
colnames(cluster3)[colnames(cluster3)=="genes"]<-"Geneid"
cluster3<-merge(cluster3,clustering_sig_genes,by='Geneid')
cluster3<-subset(cluster3,select=-c(Gene_Function.x))
```

The issue is, in her "results_table_LRT_01" object, she has 67 genes, while when I run it, I only get 33 genes. 

### Here is my code: 

```{r}
#first need to create metadata file for all hours
metadata = data.frame(sample=colnames(countsmatrix),
                condition = stringr::str_detect(pattern = ".*C.*",string = colnames(countsmatrix)),
                hour = stringr::str_replace(pattern = "(.).*",replacement="\\1",string = colnames(countsmatrix)),
 id = stringr::str_replace(pattern=".(...).*",replacement="\\1",string=colnames(countsmatrix)))

#changing TRUE and FALSE for condition to Control and Wounded
metadata$condition[str_detect(metadata$condition,"TRUE")] <- "Control"
metadata$condition[str_detect(metadata$condition,"FALSE")] <- "Wounded"

metadata <- metadata %>% column_to_rownames("sample")

#need to change condition and hour to factors
metadata$condition <- as.factor(metadata$condition)
metadata$hour<-as.factor(metadata$hour)

#LRT is the Likelihood Ratio Test, used to identify differentially expressed genes (DEGs) between control vs. treatment over time
dds_all=DESeq2::DESeqDataSetFromMatrix(countData = countsmatrix, colData = metadata, design = ~condition+hour+condition:hour)
dds_all <- estimateSizeFactors(dds_all)

#prefiltering recommended by DESeq2 Vignette (removes anything with reads below 10)
keep_all <- rowSums(counts(dds_all)) >= 10
dds_all <- dds_all[keep_all,]

#there may be an issue with Sami's code because everything is assigned as the variable "dds" which may be why she gets 67 genes instead of 33 genes for the p<0.01 LRT results? 

#sami's code:
dds_savedfile_from_old_run <- readRDS(file = "../WH_LRT_DESeq.rds")
res_savedfile_from_old_run <- results(dds_savedfile_from_old_run, alpha = 0.01)
summary(res_savedfile_from_old_run)

#DESeq2 variable for LRT test
dds_LRT<-DESeq2::DESeq(dds_all,reduced=~hour,test="LRT")

res_LRT<-results(dds_LRT)
summary(res_LRT)
# out of 21267 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 69, 0.32%
# LFC < 0 (down)     : 103, 0.48%
# outliers [1]       : 0, 0%
# low counts [2]     : 4123, 19%

#LRTs are very "generous", so a more stringent alpha value than 0.05 is needed
res_LRT_01 <- results(dds_LRT, alpha = 0.01)
summary(res_LRT_01)
# out of 21267 with nonzero total read count
# adjusted p-value < 0.01
# LFC > 0 (up)       : 11, 0.052%
# LFC < 0 (down)     : 22, 0.1%
# outliers [1]       : 0, 0%
# low counts [2]     : 3711, 17%

#annotate with Gene Function
results_LRT_table<-merge(as.data.frame(res_LRT_01),Pdam_Gene_Names_Info,by='row.names',all=FALSE)%>%
  column_to_rownames("Row.names")
dim(results_LRT_table) #21267 x 7

#filtering by a significance value of p < 0.01 for gene expression 
results_table_LRT_01 <- results_LRT_table %>% rownames_to_column("Geneid")%>% filter(padj < 0.01) #33 genes
#67 

clustering_sig_genes <- results_LRT_table %>% 
                  rownames_to_column(var = "Geneid") %>% 
                  arrange(padj) %>%
                  head(n=1000)

rld<-rlog(dds_LRT,blind=TRUE)
rld_mat <- assay(rld)
# Input is a matrix of log transformed values

rld_cor <- cor(rld_mat)    ## cor() is a base R function

head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames

#DestoryX formula for row
destroyX_row = function(es) {
  f = es
  for (row in c(1:nrow(f))){ #for each column in dataframe
    if (startsWith(rownames(f)[row], "X") == TRUE)  { #if starts with 'X' ..
      rownames(f)[row] <- substr(rownames(f)[row], 2, 100) #get rid of it
    }
  }
  assign(deparse(substitute(es)), f, inherits = TRUE) #assign corrected data to original name
}

rld_cor<-destroyX_row(rld_cor)
cluster_rlog <- rld_mat[clustering_sig_genes$Geneid, ]
#destroyX(cluster_rlog)

# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups 

clusters <- degPatterns(cluster_rlog, metadata, time="hour", col="condition")

clusters$benchmarking_curve
clusters$benchmarking
clusters$plot

#editing the plot for publication
cluster_plot<-clusters$plot+
  xlab(paste0("Hour"))+
  # theme(text=element_text(size=12))+
  # #theme(legend.key.size=unit(0.5,"cm"))+
  # theme(legend.title=element_text(size=12))+
  labs(fill="Condition")+
  scale_color_discrete(name="Condition")
  #theme_classic(base_size=12,base_family="serif")
#ggsave(filename="Cluster_Plot.png")
  
# What type of data structure is the `clusters` output?
class(clusters)

# Let's see what is stored in the `df` component
head(clusters$df)
```

I'm going to try now to run her code in her initial R script (previously I was copying and pasting sections over) to see if I can recreate her code. 
Ok I just did that and I still get 33 genes!!!! 

Now i'm thinking maybe there's something wrong with the p-adj formula or something, where I'm getting a different number of significant genes because maybe there's some different formula happening in the background or something.

I really want to try to recreate it because she got these nice 3 gene clusters, which I wanted to run through gene ontology enrichment analysis. However,
when I do it, I get nothing because there are zero genes selected from the 33 that meet whatever cutoff is default in DEGPatterns. 

```{r}
> clusters <- degPatterns(cluster_rlog, metadata, time="hour", col="condition")
#Working with 33 genes.
#Working with 0 genes after filtering: minc > 15
#Joining, by = "merge"
#Joining, by = "merge"
#Error in degPatterns(cluster_rlog, metadata, time = "hour", col = "condition") : 
  #object 'p' not found
 ```

My next thought is maybe it's the loading of the "vegan" package in mine that messes things up behind-the-scenes in R somehow, as when i do library(vegan) it has all of these notes that pop up. So I'm going to restart R and not load that package and see what happens.

Just kidding, by process of elimination all the notes that pop up (I think) are from DESeq2. 

I'll still try running the downstream code anyways and see if I still get 33 genes. I may just need to abandon ship on this LRT thing, especially since it is no longer as pertinent to the story.

Nope still get 33 genes.

Ok I'm just going to abandon ship.

### Here is the last updated code for the gene clustering analysis (I'm going to remove it from my final_pub_script.md file):

```{r}
#first need to create metadata file for all hours
metadata = data.frame(sample=colnames(countsmatrix),
                condition = stringr::str_detect(pattern = ".*C.*",string = colnames(countsmatrix)),
                hour = stringr::str_replace(pattern = "(.).*",replacement="\\1",string = colnames(countsmatrix)),
 id = stringr::str_replace(pattern=".(...).*",replacement="\\1",string=colnames(countsmatrix)))

#changing TRUE and FALSE for condition to Control and Wounded
metadata$condition[str_detect(metadata$condition,"TRUE")] <- "Control"
metadata$condition[str_detect(metadata$condition,"FALSE")] <- "Wounded"

metadata <- metadata %>% column_to_rownames("sample")

#need to change condition and hour to factors
metadata$condition <- as.factor(metadata$condition)
metadata$hour<-as.factor(metadata$hour)

#LRT is the Likelihood Ratio Test, used to identify differentially expressed genes (DEGs) between control vs. treatment over time
dds_all=DESeq2::DESeqDataSetFromMatrix(countData = countsmatrix, colData = metadata, design = ~condition+hour+condition:hour)
dds_all <- estimateSizeFactors(dds_all)

#prefiltering recommended by DESeq2 Vignette (removes anything with reads below 10)
keep_all <- rowSums(counts(dds_all)) >= 10
dds_all <- dds_all[keep_all,]

#DESeq2 variable for LRT test
dds_LRT<-DESeq2::DESeq(dds_all,reduced=~hour,test="LRT")

res_LRT<-results(dds_LRT)
summary(res_LRT)
# out of 21267 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 69, 0.32%
# LFC < 0 (down)     : 103, 0.48%
# outliers [1]       : 0, 0%
# low counts [2]     : 4123, 19%

#LRTs are very "generous", so a more stringent alpha value than 0.05 is needed
res_LRT_01 <- results(dds_LRT, alpha = 0.01)
summary(res_LRT_01)
# out of 21267 with nonzero total read count
# adjusted p-value < 0.01
# LFC > 0 (up)       : 11, 0.052%
# LFC < 0 (down)     : 22, 0.1%
# outliers [1]       : 0, 0%
# low counts [2]     : 3711, 17%

#annotate with Gene Function
results_LRT_table<-merge(as.data.frame(res_LRT_01),Pdam_Gene_Names_Info,by='row.names',all=FALSE)%>%
  column_to_rownames("Row.names")
dim(results_LRT_table) #21267 x 7

#filtering by a significance value of p < 0.01 for gene expression 
results_table_LRT_01 <- results_LRT_table %>% rownames_to_column("Geneid") %>% filter(padj < 0.01) #33 genes

clustering_sig_genes <- results_table_LRT_01 %>% 
                  arrange(padj) %>%
                  head(n=1000)

rld<-rlog(dds_LRT,blind=TRUE)
rld_mat <- assay(rld)
# Input is a matrix of log transformed values

rld_cor <- cor(rld_mat)    ## cor() is a base R function

head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames

cluster_rlog <- rld_mat[clustering_sig_genes$Geneid, ]


# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups 

clusters <- degPatterns(cluster_rlog, metadata, time="hour", col="condition")

clusters$benchmarking_curve
clusters$benchmarking
clusters$plot

#editing the plot for publication
cluster_plot<-clusters$plot+
  xlab(paste0("Hour"))+
  # theme(text=element_text(size=12))+
  # #theme(legend.key.size=unit(0.5,"cm"))+
  # theme(legend.title=element_text(size=12))+
  labs(fill="Condition")+
  scale_color_discrete(name="Condition")
  #theme_classic(base_size=12,base_family="serif")
#ggsave(filename="Cluster_Plot.png")
  
# What type of data structure is the `clusters` output?
class(clusters)

# Let's see what is stored in the `df` component
head(clusters$df)
```

Create tables for clusters
```{r}
#extract raw data from clusters
clustersdata<-as.data.frame(clusters$raw)

#cluster 1 table
cluster1<-clustersdata[clustersdata$cluster=='1',]
#clean up table and annotate with gene function
cluster1<-subset(cluster1,merge=="Control0",select=-c(merge,condition,cluster,id,hour))
#picked Control0 arbitrarily (cluster doesn't differ by condition or hour)
rownames(cluster1) <- cluster1$genes 
cluster1<-merge(cluster1,Pdam_Gene_Names_Info,by='row.names',all=TRUE)%>%
  na.omit(cluster1)
cluster1<-subset(cluster1,select=-c(Row.names,value))
#merge with cluster_sig_genes table to sort by significance
colnames(cluster1)[colnames(cluster1)=="genes"]<-"Geneid"
cluster1<-merge(cluster1,clustering_sig_genes,by='Geneid')
cluster1<-subset(cluster1,select=-c(Gene_Function.x))

#cluster 2 table
cluster2<-clustersdata[clustersdata$cluster=='2',]
#clean up table and annotate with gene function
cluster2<-subset(cluster2,merge=="Control0",select=-c(merge,condition,cluster,id,hour))
rownames(cluster2) <- cluster2$genes 
cluster2<-merge(cluster2,Pdam_Gene_Names_Info,by='row.names',all=TRUE)%>%
  na.omit(cluster2)
cluster2<-subset(cluster2,select=-c(Row.names,value))
#merge with cluster_sig_genes
colnames(cluster2)[colnames(cluster2)=="genes"]<-"Geneid"
cluster2<-merge(cluster2,clustering_sig_genes,by='Geneid')
cluster2<-subset(cluster2,select=-c(Gene_Function.x))

#cluster 3 table
cluster3<-clustersdata[clustersdata$cluster=='3',]
#clean up table and annotate with gene function
cluster3<-subset(cluster3,merge=="Control0",select=-c(merge,condition,cluster,id,hour))
rownames(cluster3) <- cluster3$genes 
cluster3<-merge(cluster3,Pdam_Gene_Names_Info,by='row.names',all=TRUE)%>%
  na.omit(cluster3)
cluster3<-subset(cluster3,select=-c(Row.names,value))
#merge with cluster_sig_genes
colnames(cluster3)[colnames(cluster3)=="genes"]<-"Geneid"
cluster3<-merge(cluster3,clustering_sig_genes,by='Geneid')
cluster3<-subset(cluster3,select=-c(Gene_Function.x))
```
