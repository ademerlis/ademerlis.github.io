---
layout: post
title: Adding gene lists to GO tables
date: '2023-05-08'
categories: [GO Analysis]
tags: [wound healing, R, Transcriptomics, GO Enrichment]
---
I want to figure out which specific genes are corresponding to each GO term -- right now it's just number of genes reported in the "Significant" column in the results table.

![image](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/baa5f1ce-d023-470e-a92f-a33a5edfe15b)

This is what I found on the avrilomics blog post: 

```{r}
> myterms = c("GO:0007610", "GO:0014070", "GO:0045910")
> mygenes <- genesInTerm(myGOdata, myterms)
> for (i in 1:length(myterms))
   {
       myterm <- myterms[i]
       mygenesforterm <- mygenes[myterm][[1]]
       mygenesforterm <- paste(mygenesforterm, collapse=',')
       print(paste("Term",myterm,"genes:",mygenesforterm))
     }
```

So here is what I've tried so far:

```{r}
GO_MF_h1_up_resultstable %>% 
  filter(Significant >= 1) %>% 
  mutate(hour=1) %>% 
  mutate(regulation_direction = "Up") %>% 
  dplyr::select(GO.ID:Fisherclassic, hour, regulation_direction) %>% 
  filter(Fisherclassic < 0.05) -> sigGO_MF_h1_up #filter for genes from this time point that are significantly associated with the GO terms (which are also significantly enriched)

length(sigGenes(GO_MF_h1_up)) #15 genes

GOterms_MF_h1_up <- sigGO_MF_h1_up %>% dplyr::select(GO.ID) %>% distinct()

GOterms_MF_h1_up <- c("GO:0042626", "GO:0015399", "GO:0016887", "GO:0003824", "GO:0022804", "GO:0005215", "GO:0140657", "GO:0008146", "GO:0016782", "GO:0017111", "GO:0016462", "GO:0016818", "GO:0016817", "GO:0022857", "GO:0008476", "GO:0016831", "GO:0008484", "GO:0016787", "GO:0016830", "GO:0051287", "GO:0016616")

genes_MF_h1_up <- genesInTerm(GO_MF_h1_up, GOterms_MF_h1_up)

for (i in 1:length(GOterms_MF_h1_up))
   {
       GOterms_MF_h1_up <- GOterms_MF_h1_up[i]
       mygenesforterm <- genes_MF_h1_up[GOterms_MF_h1_up][[1]]
       mygenesforterm <- paste(mygenesforterm, collapse=',')
       print(paste("Term",GOterms_MF_h1_up,"genes:",mygenesforterm))
     }
```

But right now I'm having trouble with these lines:
```{r}
GOterms_MF_h1_up <- sigGO_MF_h1_up %>% dplyr::select(GO.ID) %>% distinct()

genes_MF_h1_up <- genesInTerm(GO_MF_h1_up, GOterms_MF_h1_up)

for (i in 1:length(GOterms_MF_h1_up))
   {
       GOterms_MF_h1_up <- GOterms_MF_h1_up[i]
       mygenesforterm <- genes_MF_h1_up[GOterms_MF_h1_up][[1]]
       mygenesforterm <- paste(mygenesforterm, collapse=',')
       print(paste("Term",GOterms_MF_h1_up,"genes:",mygenesforterm))
     }
 ```
 I get the following error: 
 
 Error in (function (classes, fdef, mtable) :
unable to find an inherited method for function ‘genesInTerm’ for signature ‘"topGOdata", "data.frame"’

So it doesn't like that my list of GO IDs is a data frame, but I tried to manually make a list of them using a c("GO:1", "GO:2", etc) and assigning that as a variable, but then when I try to run the for loop it doesn't work at all.
