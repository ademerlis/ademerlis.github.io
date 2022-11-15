---
layout: post
title: topGO Analysis
date: '2022-11-15'
categories: Coding
tags: [Coding]
---

I am working on understanding the topGO R package (Alexa et al. 2006) to test for gene ontology enrichment from the *Pocillopora damicornis* wound healing transcriptomics study.

Resources I am currently using:
1. topGO vignette: bioconductor.org/packages/devel/bioc/vignettes/topGO/inst/doc/topGO.pdf
2. https://www.biostars.org/p/350710/
3. http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html (explains difference between algorithms and Fisher's Exact Test vs. KS test)
4. https://zhiganglu.com/post/topgo-ks-test/
5. https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html
6. https://academic.oup.com/bioinformatics/article/22/13/1600/193669 (Alexa et al. 2006 publication)
7. https://www.biostars.org/p/247636/

What I am really confused about is which algorithm and statistical test to use. I just want to make a basic bar graph showing the most significantly enriched GO terms for each hour, with the condition (control vs. wounded) as the comparison for differentially expressed genes.

It looks like from source #7 that the Fisher exact test is used when the input is gene count data, and KS test should be used if the input is p-value from DGE. I think based on the line of code assigning either a "0" or a "1" based on p-value being < or > 0.05, this would indicate that my differentially expressed genes were based on p-value, not the count data directly. So the KS test would be best. 

As for the algorithm, the only straightforward comparison I found was "classic" versus "elim", and "elim" seemed better because it took into account the GO term hierarchy (parent vs. child terms?) while the "classic" algorithm treats each GO term as an independent term (which doesn't make sense because GO terms have a hierarchy when you look at the significant node diagrams. So to not take that into account would mean that you are likely overestimating actual significance (redundancy of terms)). 

So based on my understanding of this, I will use the "elim" algorithm and the KS test. I can choose which significance level to select as the cut-off, 0.05 or 0.01. 

Also, I am not sure if I have to use all types of ontology (Biological Processes, Cellular Components, Molecular Function). The one that interests me most is BP. It looks like Mike Connelly used only BP in his publication (10.1016/j.dci.2020.103717) so that should be ok to do?

Also here is the code I am currently working on:

# topGO analysis
```{r}
#following Michael Connelly's code ([github](https://github.com/michaeltconnelly/EAPSI_Pocillopora_LPS/blob/master/Rmd/LPS_topGO_DESeq2_Pdam.Rmd))

#Input required *P. damicornis* GO annotation data and construct Gene-to-GO object for custom annotation mapping
GO_geneID<-readMappings(file="../2022 updates with Sami Beasley/pdam_genome_GOgenes.txt") #Mike made this file

geneID_GO <- inverseList(GO_geneID)
str(head(geneID_GO))

#Generate gene universe and GO universe from Gene-to-GO and GO-to-Gene objects
geneNames <- names(geneID_GO)
str(head(geneNames))

GONames <- names(GO_geneID)
str(head(GONames))
```

#hour 0 no annotations so no GO analysis

## topGO analysis Hour 1 Biological Processes (BP)
```{r}
#label differentially expressed genes as '1' and non-significant genes as '0'
res_h1_GO<-as.data.frame(res_h1)
res_h1_GO<- res_h1_GO%>%mutate(Group=case_when(padj<0.05~"1",padj>=0.05~"0"))
res_h1_GO<-na.omit(res_h1_GO)

#create the gene universe
gene_universe_h1<-as.numeric(res_h1_GO$Group)
gene_universe_h1<-factor(gene_universe_h1)
names(gene_universe_h1) <- rownames(res_h1_GO) 

#create topGO data object
GO_BP_h1 <-new("topGOdata",
                ontology="BP",
                allGenes=gene_universe_h1,
                nodeSize=10,
                annotationFun=annFUN.gene2GO,
                gene2GO = geneID_GO)

#two algorithms: classic (each GO term is tested independently, not taking the GO hierarchy into account) or elim (this method processes the GO terms by traversing the GO hierarchy from bottom to top, ie. it first assesses the most specific (bottom-most) GO terms, and proceeds later to more general (higher) GO terms. When it assesses a higher (more general) GO term, it discards any genes that are annotated with significantly enriched descendant GO terms (considered significant using a pre-defined P-value threshold). This method does tend to miss some true positives at higher (more general) levels of the GO hierarchy.)
#Source: http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html 

#Fisher's exact test to identify enriched gene ontology terms
GO_BP_h1_resultFisher<-runTest(GO_BP_h1,algorithm = "elim", statistic = "fisher")

#generate results table
GO_BP_h1_resultstable<-GenTable(GO_BP_h1, 
                   elimFisher = GO_BP_h1_resultFisher,
                   topNodes = 50)

showSigOfNodes(GO_BP_h1,score(GO_BP_h1_resultFisher),firstSigNodes = 5,useInfo = 'all')

GO_BP_h1_resultstable %>% filter(elimFisher < 0.01)
```
