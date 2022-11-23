---
layout: post
title: GO vs KEGG Analysis
date: '2022-11-22'
categories: Coding
tags: [Coding]
---

In doing the gene ontology enrichment analysis for biological processes (comparing control vs. wounded at each time point), we don't get any meaningful GO terms that pop up.

```{r}
#following Michael Connelly's code (github)

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

Hour 0 had no annotated significant genes so we skip that (GO analysis isn't useful for unknown function genes)

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

```{r}
#label differentially expressed genes as '1' and non-significant genes as '0'
res_h2_GO<-as.data.frame(res_h2)
res_h2_GO<- res_h2_GO%>%mutate(Group=case_when(padj<0.05~"1",padj>=0.05~"0"))
res_h2_GO<-na.omit(res_h2_GO)

#create the gene universe
gene_universe_h2<-as.numeric(res_h2_GO$Group)
gene_universe_h2<-factor(gene_universe_h2)
names(gene_universe_h2) <- rownames(res_h2_GO) 

#create topGO data object
GO_BP_h2 <-new("topGOdata",
                ontology="BP",
                allGenes=gene_universe_h2,
                nodeSize=10,
                annotationFun=annFUN.gene2GO,
                gene2GO = geneID_GO)

# tests to identify enriched gene ontology terms
GO_BP_h2_resultFisher01<-runTest(GO_BP_h2,algorithm = "weight01", statistic = "fisher")
GO_BP_h2_resultFisherclassic<-runTest(GO_BP_h2,algorithm = "classic", statistic = "fisher")
GO_BP_h2_resultclassicKS <- runTest(GO_BP_h2, algorithm = "classic", statistic = "ks")
GO_BP_h2_resultelimKS <- runTest(GO_BP_h2, algorithm = "elim", statistic = "ks")


#generate results table
GO_BP_h2_resultstable<-GenTable(GO_BP_h2, 
                                Fisher01 = GO_BP_h2_resultFisher01,
                                Fisherclassic = GO_BP_h2_resultFisherclassic,
                                classicKS = GO_BP_h2_resultclassicKS,
                                elimKS = GO_BP_h2_resultelimKS,
                                topNodes = 50)

showSigOfNodes(GO_BP_h2,score(GO_BP_h2_resultelimKS),firstSigNodes = 5,useInfo = 'all')

GO_BP_h2_resultstable %>% filter(Significant >= 1) %>% mutate(hour=2)-> sigGO_BP_h2 #filter for genes from this time point that are significantly associated with the GO terms (which are also significantly enriched)
```

```{r}
#label differentially expressed genes as '1' and non-significant genes as '0'
res_h4_GO<-as.data.frame(res_h4)
res_h4_GO<- res_h4_GO%>%mutate(Group=case_when(padj<0.05~"1",padj>=0.05~"0"))
res_h4_GO<-na.omit(res_h4_GO)

#create the gene universe
gene_universe_h4<-as.numeric(res_h4_GO$Group)
gene_universe_h4<-factor(gene_universe_h4)
names(gene_universe_h4) <- rownames(res_h4_GO) 

#create topGO data object
GO_BP_h4 <-new("topGOdata",
                ontology="BP",
                allGenes=gene_universe_h4,
                nodeSize=10,
                annotationFun=annFUN.gene2GO,
                gene2GO = geneID_GO)

# tests to identify enriched gene ontology terms
GO_BP_h4_resultFisher01<-runTest(GO_BP_h4,algorithm = "weight01", statistic = "fisher")
GO_BP_h4_resultFisherclassic<-runTest(GO_BP_h4,algorithm = "classic", statistic = "fisher")
GO_BP_h4_resultclassicKS <- runTest(GO_BP_h4, algorithm = "classic", statistic = "ks")
GO_BP_h4_resultelimKS <- runTest(GO_BP_h4, algorithm = "elim", statistic = "ks")


#generate results table
GO_BP_h4_resultstable<-GenTable(GO_BP_h4, 
                                Fisher01 = GO_BP_h4_resultFisher01,
                                Fisherclassic = GO_BP_h4_resultFisherclassic,
                                classicKS = GO_BP_h4_resultclassicKS,
                                elimKS = GO_BP_h4_resultelimKS,
                                topNodes = 50)

showSigOfNodes(GO_BP_h4,score(GO_BP_h4_resultelimKS),firstSigNodes = 5,useInfo = 'all')

GO_BP_h4_resultstable %>% filter(Significant >= 1) %>% mutate(hour=4)-> sigGO_BP_h4 #filter for genes from this time point that are significantly associated with the GO terms (which are also significantly enriched)
```

```{r}
full_join(sigGO_BP_h1, sigGO_BP_h2) %>% full_join(., sigGO_BP_h4) %>%  #109 GO terms total
  filter(Fisher01 < 0.05) #11 GO terms now

full_join(sigGO_BP_h1, sigGO_BP_h2) %>% full_join(., sigGO_BP_h4) %>%  #109 GO terms total
  filter(Fisherclassic < 0.05) #16 GO terms now

full_join(sigGO_BP_h1, sigGO_BP_h2) %>% full_join(., sigGO_BP_h4) %>%
  filter(Fisherclassic < 0.05) %>% 
  mutate(GO_term=str_c(GO.ID,Term, sep=" ")) %>% 
  mutate(Fisherclassic=as.numeric(Fisherclassic), hour=as.factor(hour)) %>% 
  ggplot(., aes(x=reorder(GO_term, -log10(Fisherclassic)), y=-log10(Fisherclassic), fill = hour)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  theme_classic() +
  coord_flip() +
  labs(y="-log10(p-value)", x="Biological Process") + 
  facet_wrap(~hour)
#ggsave("significant_GO_BP_perhour.pdf")
```

<img width="695" alt="Screen Shot 2022-11-23 at 9 17 14 AM" src="https://user-images.githubusercontent.com/56000927/203569199-13c77704-1b99-4623-a77a-6e3c95235cfc.png">


We're going to try this now for molecular function (when Sami did it, she got terms related to peroxidase activity, which is immune related).

```{r}
#create topGO data object
GO_MF_h1 <-new("topGOdata",
                ontology="MF",
                allGenes=gene_universe_h1,
                nodeSize=10,
                annotationFun=annFUN.gene2GO,
                gene2GO = geneID_GO)

# tests to identify enriched gene ontology terms
GO_MF_h1_resultFisher01<-runTest(GO_MF_h1,algorithm = "weight01", statistic = "fisher")
GO_MF_h1_resultFisherclassic<-runTest(GO_MF_h1,algorithm = "classic", statistic = "fisher")
GO_MF_h1_resultclassicKS <- runTest(GO_MF_h1, algorithm = "classic", statistic = "ks")
GO_MF_h1_resultelimKS <- runTest(GO_MF_h1, algorithm = "elim", statistic = "ks")


#generate results table
GO_MF_h1_resultstable<-GenTable(GO_MF_h1, 
                                Fisher01 = GO_MF_h1_resultFisher01,
                                Fisherclassic = GO_MF_h1_resultFisherclassic,
                                classicKS = GO_MF_h1_resultclassicKS,
                                elimKS = GO_MF_h1_resultelimKS,
                                topNodes = 50)

showSigOfNodes(GO_MF_h1,score(GO_MF_h1_resultelimKS),firstSigNodes = 5,useInfo = 'all')

GO_MF_h1_resultstable %>% filter(Significant >= 1) %>% mutate(hour=1)-> sigGO_MF_h1 #filter for genes from this time point that are significantly associated with the GO terms (which are also significantly enriched)

#create topGO data object
GO_MF_h2 <-new("topGOdata",
                ontology="MF",
                allGenes=gene_universe_h2,
                nodeSize=10,
                annotationFun=annFUN.gene2GO,
                gene2GO = geneID_GO)

# tests to identify enriched gene ontology terms
GO_MF_h2_resultFisher01<-runTest(GO_MF_h2,algorithm = "weight01", statistic = "fisher")
GO_MF_h2_resultFisherclassic<-runTest(GO_MF_h2,algorithm = "classic", statistic = "fisher")
GO_MF_h2_resultclassicKS <- runTest(GO_MF_h2, algorithm = "classic", statistic = "ks")
GO_MF_h2_resultelimKS <- runTest(GO_MF_h2, algorithm = "elim", statistic = "ks")


#generate results table
GO_MF_h2_resultstable<-GenTable(GO_MF_h2, 
                                Fisher01 = GO_MF_h2_resultFisher01,
                                Fisherclassic = GO_MF_h2_resultFisherclassic,
                                classicKS = GO_MF_h2_resultclassicKS,
                                elimKS = GO_MF_h2_resultelimKS,
                                topNodes = 50)

showSigOfNodes(GO_MF_h2,score(GO_MF_h2_resultelimKS),firstSigNodes = 5,useInfo = 'all')

GO_MF_h2_resultstable %>% filter(Significant >= 1) %>% mutate(hour=2)-> sigGO_MF_h2 #filter for genes from this time point that are significantly associated with the GO terms (which are also significantly enriched)

#create topGO data object
GO_MF_h4 <-new("topGOdata",
                ontology="MF",
                allGenes=gene_universe_h4,
                nodeSize=10,
                annotationFun=annFUN.gene2GO,
                gene2GO = geneID_GO)

# tests to identify enriched gene ontology terms
GO_MF_h4_resultFisher01<-runTest(GO_MF_h4,algorithm = "weight01", statistic = "fisher")
GO_MF_h4_resultFisherclassic<-runTest(GO_MF_h4,algorithm = "classic", statistic = "fisher")
GO_MF_h4_resultclassicKS <- runTest(GO_MF_h4, algorithm = "classic", statistic = "ks")
GO_MF_h4_resultelimKS <- runTest(GO_MF_h4, algorithm = "elim", statistic = "ks")


#generate results table
GO_MF_h4_resultstable<-GenTable(GO_MF_h4, 
                                Fisher01 = GO_MF_h4_resultFisher01,
                                Fisherclassic = GO_MF_h4_resultFisherclassic,
                                classicKS = GO_MF_h4_resultclassicKS,
                                elimKS = GO_MF_h4_resultelimKS,
                                topNodes = 50)

showSigOfNodes(GO_MF_h4,score(GO_MF_h4_resultelimKS),firstSigNodes = 5,useInfo = 'all')

GO_MF_h4_resultstable %>% filter(Significant >= 1) %>% mutate(hour=4)-> sigGO_MF_h4 #filter for genes from this time point that are significantly associated with the GO terms (which are also significantly enriched)

full_join(sigGO_MF_h1, sigGO_MF_h2) %>% full_join(., sigGO_MF_h4) %>%  #144 GO terms total
  filter(Fisher01 < 0.05) #14 GO terms now

full_join(sigGO_MF_h1, sigGO_MF_h2) %>% full_join(., sigGO_MF_h4) %>%  #144 GO terms total
  filter(Fisherclassic < 0.05) #25 GO terms now

full_join(sigGO_MF_h1, sigGO_MF_h2) %>% full_join(., sigGO_MF_h4) %>%
  filter(Fisherclassic < 0.05) %>% 
  mutate(GO_term=str_c(GO.ID,Term, sep=" ")) %>% 
  mutate(Fisherclassic=as.numeric(Fisherclassic), hour=as.factor(hour)) %>% 
  ggplot(., aes(x=reorder(GO_term, -log10(Fisherclassic)), y=-log10(Fisherclassic), fill = hour)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  theme_classic() +
  coord_flip() +
  labs(y="-log10(p-value)", x="Molecular Function") + 
  facet_wrap(~hour)
#ggsave("significant_GO_MF_perhour.pdf")
```

<img width="688" alt="Screen Shot 2022-11-23 at 9 36 17 AM" src="https://user-images.githubusercontent.com/56000927/203573638-cc133ac0-b648-43c9-982e-d77403a3e5cc.png">

We get some cool terms related to immune system activity!!! 


