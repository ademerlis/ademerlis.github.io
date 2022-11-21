---
layout: post
title: PCA plot significance
date: '2022-11-21'
categories: Coding
tags: [Coding]
---

I have made my PCA plots for each hour of the differential expression analysis (Figure below).

![PCA_woundhealing](https://user-images.githubusercontent.com/56000927/203095584-e9d04fa2-eb9d-484d-9b73-70773906744c.jpg)

I want to add statistical tests and significance to support the differential gene expression at each time point (suggestion from Dr. Kevin Wong - thank you!).

Kevin sent me his Physiology Analysis Rmd for a paper he's currently working on, which has this PCA in it with PERMANOVAs included:

<img width="718" alt="Screen Shot 2022-11-21 at 10 49 21 AM" src="https://user-images.githubusercontent.com/56000927/203098733-358844a0-bde7-44a7-bbb2-0bc8456178ac.png">

I want to also do a PERMANOVA for my PCA plots.

What I've got so far:
Example of DDS object creation for hour 0 (control vs. wounded):

```{r}
countdata_0 <- countsmatrix %>% dplyr::select(`0301C1`:`0306Z`)

metadata_0 = data.frame(sample=colnames(countdata_0),
                condition = stringr::str_detect(pattern = ".*C.*",string = colnames(countdata_0)),
                hour = stringr::str_replace(pattern = "(.).*",replacement="\\1",string = colnames(countdata_0)),
 id = stringr::str_replace(pattern="(...).*",replacement="\\1",string=colnames(countdata_0)))

#changing TRUE and FALSE for condition to Control and Wounded
metadata_0$condition[str_detect(metadata_0$condition,"TRUE")] <- "Control"
metadata_0$condition[str_detect(metadata_0$condition,"FALSE")] <- "Wounded"

metadata_0 <- metadata_0 %>% column_to_rownames("sample")

metadata_0$condition <- as.factor(metadata_0$condition)

#DESeq2 from a counts matrix; requires count data, metadata, and the experimental design
dds_0=DESeq2::DESeqDataSetFromMatrix(countData = countdata_0, colData = metadata_0, design = ~condition)
dds_0 <- estimateSizeFactors(dds_0)

#prefiltering recommended by DESeq2 Vignette (removes anything with reads below 10)
keep_0 <- rowSums(counts(dds_0)) >= 10
dds_0 <- dds_0[keep_0,]

dds_0<-DESeq(dds_0)

res_h0<-results(dds_0)
summary(res_h0,alpha=0.05)
#5 downregulated genes (3 outliers), no upregulated genes

#identifying significant genes
resSig_h0<-res_h0[which(res_h0$padj<0.05), ]

#annotate results to include Gene Function
anno_h0<-merge(as.data.frame(resSig_h0),Pdam_Gene_Names_Info,by='row.names',all=TRUE)
anno_h0<-na.omit(anno_h0)
```

# Principal component analysis plot and Scree plots
```{r}
#need to transform the data for plotting
dds_vst0<- vst(dds_0,blind=FALSE)
plotPCA(dds_vst0)
#plotPCA is one function to do it, or you can use the "prcomp" function, which lets you create a scree plot
pca_h0 <- prcomp(t(assay(dds_vst0)))
fviz_eig(pca_h0)
```

## PCA figures for manuscript
```{r}
#plotting the PCA in ggplot
pca12_0 <- plotPCA(dds_vst0,intgroup=c("condition"),returnData = TRUE)
pca_0 = ggplot(pca12_0, aes(PC1,PC2,shape=condition,color=condition)) + 
  geom_point(size=3) +  
  xlab(paste0("PC1 (49%)")) +
  ylab(paste0("PC2 (24%)")) +
  theme(legend.position="right")  + 
  theme(text = element_text(size=12))  + 
  theme(legend.key.size = unit(0.5, "cm")) + 
  theme(legend.title=element_text(size=12)) +
  scale_shape_discrete(name="Condition") +
  scale_color_discrete(name="Condition") +
  theme_classic(base_size=12,base_family="serif")
  ```
  
  Now I'm having trouble with the PERMANOVA part.
  
 # This is the section of code Kevin used: 
  ```{r}
  master <- join_all(dfs, by=c("Fragment.ID", "Day", "Group"))
  master.pca <- master
  coral_info<-master.pca[c(2,3)] #columns 2 and 3 are "Day" and "Group" in his metadata file

#Examine PERMANOVA results.  
# scale data
vegan <- scale(master.pca[c(4:18)]) #these columns he created as calculations for physiology metrics (i.e. ChlC2.ugcell,  Carb.mgcell, Protein.mgcell)

# PerMANOVA 
permanova<-adonis(vegan ~ Day*Group, data = master.pca, method='eu')
z_pca<-permanova$aov.tab
z_pca

## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## Day        2    174.38  87.189  9.0242 0.26421  0.001 ***
## Group      2     78.41  39.205  4.0577 0.11880  0.002 ** 
## Day:Group  4     59.39  14.848  1.5368 0.08999  0.069 .  
## Residuals 36    347.82   9.662         0.52700           
## Total     44    660.00                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
 ```
