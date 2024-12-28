---
layout: post
title: topGO Analysis
date: '2022-11-15'
categories: [Statistics for GO Analysis]
tags: [Transcriptomics, wound healing, Pocillopora damicornis, topGO, R, Statistics, GO Enrichment]
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

Ok, after following my own train of thought from above (using elim + KS), the results from each hour were sparse and very general GO terms. In looking at the resource links a bit more, one said that the selection of the test and algorithm that you end up using depends on how "useful" the GO terms are that it gives you in the results. So it sounds like you can pick whatever combination works best for you.

So in the code below, I run 4 different combinations of algorithms and statistical tests (classic+Fisher, weight01+Fisher, elim+KS, classic+KS) to see what comes out as significant and which GO terms are interesting.

Here is an example just from hour 1. 

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

# tests to identify enriched gene ontology terms
GO_BP_h1_resultFisher01<-runTest(GO_BP_h1,algorithm = "weight01", statistic = "fisher")
GO_BP_h1_resultFisherclassic<-runTest(GO_BP_h1,algorithm = "classic", statistic = "fisher")
GO_BP_h1_resultclassicKS <- runTest(GO_BP_h1, algorithm = "classic", statistic = "ks")
GO_BP_h1_resultelimKS <- runTest(GO_BP_h1, algorithm = "elim", statistic = "ks")


#generate results table
GO_BP_h1_resultstable<-GenTable(GO_BP_h1, 
                                Fisher01 = GO_BP_h1_resultFisher01,
                                Fisherclassic = GO_BP_h1_resultFisherclassic,
                                classicKS = GO_BP_h1_resultclassicKS,
                                elimKS = GO_BP_h1_resultelimKS,
                                topNodes = 50)

showSigOfNodes(GO_BP_h1,score(GO_BP_h1_resultelimKS),firstSigNodes = 5,useInfo = 'all')

GO_BP_h1_resultstable %>% filter(Significant >= 1) %>% mutate(hour=1)-> sigGO_BP_h1 #filter for genes from this time point that are significantly associated with the GO terms (which are also significantly enriched)
```
So the results go from the first screenshot (just the elim+KS test) to the second one (multiple test combinations
![Screen Shot 2022-11-15 at 11 07 00 AM](https://user-images.githubusercontent.com/56000927/201968433-89a7063e-f56a-468b-9afa-33a0e0dd7e9f.png)
![Screen Shot 2022-11-15 at 11 07 19 AM](https://user-images.githubusercontent.com/56000927/201968471-2c7c7317-270c-4678-ac11-b9ba680f5106.png)

And as you can see, the second one also has more specific GO terms, which I think is due to the use of the Fisher exact test, which is not accounting for GO term hierarchy. So maybe actually that is useful to treat GO terms as independent, otherwise when using the KS test, it only picks the top parent GO term as the signficant one, which isn't very descriptive (could be "cellular process" or something). 

Ok so based on all this I think I will go forward with Fisher exact test and either the classic or the weight01 algorithm.
