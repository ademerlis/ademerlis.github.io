---
layout: post
title: PCA plot significance
date: '2022-11-21'
categories: [Trialing code]
tags: [R, PCA, wound healing, DESeq2, PERMANOVA]
---

I have made my PCA plots for each hour of the differential expression analysis (Figure below).

![PCA_woundhealing](https://user-images.githubusercontent.com/56000927/203095584-e9d04fa2-eb9d-484d-9b73-70773906744c.jpg)

I want to add statistical tests and significance to support the differential gene expression at each time point (suggestion from Dr. Kevin Wong - thank you!).

Kevin sent me his Physiology Analysis Rmd for a paper he's currently working on, which has this PCA in it with PERMANOVAs included:

<img width="718" alt="Screen Shot 2022-11-21 at 10 49 21 AM" src="https://user-images.githubusercontent.com/56000927/203098733-358844a0-bde7-44a7-bbb2-0bc8456178ac.png">

I want to also do a PERMANOVA for my PCA plots.

What I've got so far:

## Example of DDS object creation for hour 0 (control vs. wounded):

```{r}
countdata_0 <- countsmatrix %>% dplyr::select(`0301C1`:`0306Z`)

metadata_0 = data.frame(sample=colnames(countdata_0),
                condition = stringr::str_detect(pattern = ".*C.*",string = colnames(countdata_0)),
                hour = stringr::str_replace(pattern = "(.).*",replacement="\\1",string = colnames(countdata_0)),
 id = stringr::str_replace(pattern=".(...).*",replacement="\\1",string=colnames(countdata_0)))

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

## Principal component analysis plot and Scree plots
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
  
 ## This is the section of code Kevin used: 
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
 
 I think I need to make a dataset like this
 
 | Sample ID  | Condition | Pdam0000001 |
| ------------- | ------------- | ------------- |
| 301 | Control  | 6 |
| 302 | Control  | 18|
| 303 | Control  | 12 |
| 304 | Wounded  | 12 |
| 305 | Wounded  | 20 |
| 306 | Wounded  | 43 |

Ok, when I rearranged the data to look like that, the permanova didn't work.

# PERMANOVA 
```{r}
#first need to restructure results matrix so that it includes each row as a sample with the condition, and then each column is a gene
hour0_countsmatrix <- assay(dds_0)
t(hour0_countsmatrix) -> hour0_countsmatrix 
#we use the dds object instead of the original countdata_0 matrix because this has filtered out genes with low counts (<10)
PCA.h0.countsdata <- merge(metadata_0, hour0_countsmatrix, by='row.names', all=TRUE)
dim(PCA.h0.countsdata)
# scale data
vegan <- scale(PCA.h0.countsdata[c(5:19451)]) #we just want to scale the gene counts

# PERMANOVA 
permanova<-adonis2(vegan ~ condition, data = PCA.h0.countsdata, method='eu', na.rm = TRUE, nperm = 999)

```
<img width="424" alt="Screen Shot 2022-11-21 at 1 13 18 PM" src="https://user-images.githubusercontent.com/56000927/203129685-78020370-12f8-463d-8ad2-8cf864be8276.png">

So I think maybe it has to be formatted a different way. 

 | Sample ID  | Condition | Gene | Count |
| ------------- | ------------- | ------------- | ------------- |
| 301 | Control  | Pdam0001 | 6 |
| 302 | Control  | Pdam0001 |18|
| 303 | Control  | Pdam0001 |12 |
| 304 | Wounded  | Pdam0001 |12 |
| 305 | Wounded  | Pdam0001 |20 |
| 306 | Wounded  | Pdam0001 |43 |

And then the design for the PERMANOVA would actually be count ~ gene*condition

Let's try that:

```{r}
#first need to restructure results matrix so that it includes each row as a sample with the condition, and then each column is a gene
#we use the dds object instead of the original countdata_0 matrix because this has filtered out genes with low counts (<10)
hour0_countsmatrix <- assay(dds_0)
head(hour0_countsmatrix) #first gene is pdam_00021773
tail(hour0_countsmatrix) #last gene is pdam_00025493
t(hour0_countsmatrix) %>% as.data.frame() -> hour0_countsmatrix 

PCA.h0.countsdata <- merge(metadata_0, hour0_countsmatrix, by='row.names', all=TRUE)
PCA.h0.countsdata %>% 
  pivot_longer(cols = pdam_00021773:pdam_00025493, names_to = "Gene ID", values_to = "Count") %>% as.matrix()->PCA.h0.countsdata_longformat

# PERMANOVA 
permanova<-adonis2(Count ~ `Gene ID`*condition, data = PCA.h0.countsdata_longformat, method='eu')
```
Doesn't work.
Error:Error in eval(YVAR, environment(formula), globalenv()) : 
  object 'Count' not found

This is what the long format data frame looks like:

<img width="514" alt="Screen Shot 2022-11-21 at 1 42 30 PM" src="https://user-images.githubusercontent.com/56000927/203134673-d52ae072-43c8-4fcd-94d9-65dcd8b91704.png">

I'm so confused. I thought maybe i would want Gene and Condition and Gene:Condition as the terms for the PERMANOVA. I need to look up adonis2 and what it accepts in the formula.

From [adonis CRAN file](https://search.r-project.org/CRAN/refmans/vegan/html/adonis.html)
"The left-hand side (LHS) of the formula must be either a community data matrix or a dissimilarity matrix, e.g., from vegdist or dist."

So that's why Kevin used "vegan" as his input variable in the formula, which is a matrix of just his dependent variables.

Adonis/Vegan is essentially running some sort of dissimilarity assessment, so it needs a matrix to run it on. So you can't put all the counts into one column.

I think I'm back to square one, where I actually got a result (133-150). I think it worked but it didn't give me significant results. I ran it for each hour as well, here is that code and screenshots of the results below.

# PERMANOVA for hour 0
```{r}
#we use the dds object instead of the original countdata_0 matrix because this has filtered out genes with low counts (<10)
hour0_countsmatrix <- assay(dds_0)
t(hour0_countsmatrix) -> hour0_countsmatrix 
PCA.h0.countsdata <- merge(metadata_0, hour0_countsmatrix, by='row.names', all=TRUE)

# scale data
vegan <- scale(PCA.h0.countsdata[c(5:19451)]) #we just want to scale the gene counts

# PERMANOVA 
permanova<-adonis2(vegan ~ condition, data = PCA.h0.countsdata, method='eu', na.rm=TRUE, nperm = 999)
permanova
```
<img width="768" alt="Screen Shot 2022-11-21 at 3 19 19 PM" src="https://user-images.githubusercontent.com/56000927/203150692-50d84201-165b-46cc-946f-cdf85f904a1b.png">

# PERMANOVA for hour 1
```{r}
#we use the dds object instead of the original countdata_0 matrix because this has filtered out genes with low counts (<10)
hour1_countsmatrix <- assay(dds_1)
t(hour1_countsmatrix) -> hour1_countsmatrix 
PCA.h1.countsdata <- merge(metadata_1, hour1_countsmatrix, by='row.names', all=TRUE)

# scale data
vegan <- scale(PCA.h1.countsdata[c(5:19536)]) #we just want to scale the gene counts

# PERMANOVA 
permanova<-adonis2(vegan ~ condition, data = PCA.h1.countsdata, method='eu', na.rm=TRUE, nperm = 999)
permanova
```

<img width="782" alt="Screen Shot 2022-11-21 at 3 19 37 PM" src="https://user-images.githubusercontent.com/56000927/203150747-0e2676bd-1c06-443b-9578-c4492268519b.png">


# PERMANOVA for hour 2
```{r}
#we use the dds object instead of the original countdata_0 matrix because this has filtered out genes with low counts (<10)
hour2_countsmatrix <- assay(dds_2)
t(hour2_countsmatrix) -> hour2_countsmatrix 
PCA.h2.countsdata <- merge(metadata_2, hour2_countsmatrix, by='row.names', all=TRUE)

# scale data
vegan <- scale(PCA.h2.countsdata[c(5:19779)]) #we just want to scale the gene counts

# PERMANOVA 
permanova<-adonis2(vegan ~ condition, data = PCA.h2.countsdata, method='eu', na.rm=TRUE, nperm = 999)
permanova
```
<img width="766" alt="Screen Shot 2022-11-21 at 3 19 59 PM" src="https://user-images.githubusercontent.com/56000927/203150802-26b2f259-3f6f-4838-be0a-4c4f97e3b762.png">

# PERMANOVA for hour 4
```{r}
#we use the dds object instead of the original countdata_0 matrix because this has filtered out genes with low counts (<10)
hour4_countsmatrix <- assay(dds_4)
t(hour4_countsmatrix) -> hour4_countsmatrix 
PCA.h4.countsdata <- merge(metadata_4, hour4_countsmatrix, by='row.names', all=TRUE)

# scale data
vegan <- scale(PCA.h4.countsdata[c(5:19773)]) #we just want to scale the gene counts

# PERMANOVA 
permanova<-adonis2(vegan ~ condition, data = PCA.h4.countsdata, method='eu', na.rm=TRUE, nperm = 999)
permanova
```

<img width="787" alt="Screen Shot 2022-11-21 at 3 23 15 PM" src="https://user-images.githubusercontent.com/56000927/203151271-c5ebae98-4616-4f2a-a78a-591f970910d7.png">

I also tried running it without the scale() function for hour 4 and I got a less significant result.
<img width="766" alt="Screen Shot 2022-11-21 at 3 20 45 PM" src="https://user-images.githubusercontent.com/56000927/203150926-309ed260-85ba-4b61-93d6-72e8e74b9093.png">

I'm going to abandon this analysis for the time being because I don't think this is super necessary for the PCA plots anyways.
