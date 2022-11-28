---
layout: post
title: GO Analysis
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

Is there a way to make this figure look better...

I think I want to try to make a dotplot, like in [Traylor-Knowles et al. 2021](https://www.frontiersin.org/articles/10.3389/fmars.2021.681563/full)

<img width="900" alt="Screen Shot 2022-11-23 at 9 42 20 AM" src="https://user-images.githubusercontent.com/56000927/203574995-189a090c-79ce-4fe7-ae3d-25ee4f81c22f.png">

I followed Mel's code from the 2021 paper [on Github](https://github.com/ademerlis/sctld_transcriptomics_2021/blob/main/SCTLD_BiNGO_Analysis_May2021.Rmd)

Here is my code:
```{r}
#dot plot
full_join(sigGO_MF_h1, sigGO_MF_h2) %>% full_join(., sigGO_MF_h4) %>%
  filter(Fisherclassic < 0.05) %>% 
  mutate(Fisherclassic = as.numeric(Fisherclassic)) %>% 
ggplot(data=., aes(x=-log(Fisherclassic), y=reorder(Term, -Fisherclassic))) +
    geom_point(aes(size=Annotated, color=-log(Fisherclassic))) +
    scale_color_viridis(option="plasma") +
  labs(x="-log P-Value", y="Gene Ontology Term Description", title=paste("Significant Molecular Function"), color = "-log(P-value)") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +
  theme_bw()
```

<img width="694" alt="Screen Shot 2022-11-23 at 11 15 01 AM" src="https://user-images.githubusercontent.com/56000927/203595507-76663c09-9243-43a5-8afb-6dca924b8a10.png">

What's annoying me is that some of the GO terms are cut-off, and no matter how I try to scale or wrap the text on the Y axis, nothing helps. This is because I realized it's actually the annotation itself from the topGO function. It must have something where after a certain number of characters it just automatically adds "..." to the end. I have been googling and I can't figure out how to extract the whole term. Unfortunately I think I need to manually edit it in Illustrator. I can google the GO terms and get the same result using [QuickGO](https://www.ebi.ac.uk/QuickGO/search/GO:0008476) and just searching the term and picking the most annotated one. They seem to all match so far.

I also don't like how the size for the annotation circles is like super high (1000-3000) but it seems like a lot of the values are small. I am going to try to fix the range of circles given so that more sizes are provided.

I couldn't figure out how to add more circles to the Annotated key, so I just changed the range of sizes on the points on the plot themselves so they were scaled a bit larger overall. 

I also realized I needed to separate it by hour, so I added that as well.

```{r}
#dot plot
full_join(sigGO_MF_h1, sigGO_MF_h2) %>% full_join(., sigGO_MF_h4) %>%
  filter(Fisherclassic < 0.05) %>% 
  mutate(Fisherclassic = as.numeric(Fisherclassic)) %>% 
  mutate(Annotated = as.numeric(Annotated)) %>% 
ggplot(data=., aes(x=-log(Fisherclassic), y=reorder(Term, -Fisherclassic))) +
    geom_point(aes(size=Annotated, color=-log(Fisherclassic))) +
    scale_color_viridis(option="plasma") +
  labs(x="-log P-Value", y="Gene Ontology Term Description", title=paste("Significant Molecular Functions"), color = "-log(P-value)") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
  theme_bw() +
  facet_grid(vars(hour), scales = "free") + 
  scale_size(range = c(2,6))
ggsave("significantMF_GO_perhour.pdf")
```

<img width="694" alt="Screen Shot 2022-11-23 at 11 36 18 AM" src="https://user-images.githubusercontent.com/56000927/203600181-4ffed145-dff4-4d93-904d-5afc72a36015.png">

I need to bring this into Illustrator now and edit it manually. 
