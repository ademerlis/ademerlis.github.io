---
layout: post
title: Re-running TopGO with Up vs. Downregulated genes split up
date: '2023-02-08'
categories: [GO Analysis]
tags: [topGO, R, wound healing, Transcriptomics, Fisher Exact Test, GO Enrichment]
---

I just checked and I did use Fisher's Exact Test (re: last post), so I don't need to redo that. But I need to redo a couple things. First, I need to split up the Up and Downregulated gene lists and run TopGO separately. Next, I need to make sure that i filtered L2FC < or > 0, instead of 2 which I did for the volcano plots. 2 is very, very strict, and all the papers I'm reading count differential expression as p-adjusted < 0.05 and L2FC < or > 0. 
My annotated gene lists that I saved from each DESeq2 object for each hour are filtered based on p-adjusted < 0.05, but no specific L2FC. So that's good, that means the only thing I need to change is the volcano plots. Although I did a lot of manual editing in Illustrator, so maybe I can just leave that as is.
I just need to make sure that for TopGO, I count any upregulated genes as L2FC >0, and downregulated genes as L2FC < 0.

Ok, so in looking at Mike Connelly's paper (Connelly et al. 2020 - LPS and Pdam paper), he has a GO table that he published where he includes the following: 

![Screen Shot 2023-05-03 at 11 20 56 AM](https://user-images.githubusercontent.com/56000927/235961836-50e1c614-9b82-467e-bda0-5d688af02ffe.png)

In the methods he says: "Gene ontology (GO) enrichment was completed using topGO (Alexa et al., 2006). GO terms corresponding to the “Biological Process” (BP) were annotated to genes in the P. damicornis genome and tested for enrichment using the “weight01” algorithm in topGO, with significant enrichment determined using Fisher's exact tests and a significance threshold of p < 0.05."

So that far-right column that says "p-value" corresponds to the significant p-values he filtered out from the "weight01" test. 

I decided to use the "classic" algorithm for mine, so I'll filter for significance there and create a table with that.

I also think that the column for "annotated" genes is not useful to report, because from my understanding, it is the number of genes from the gene universe you inputted into TopGO that align to a given GO term. But all we really care about are the significant DGEs and whether they have a significant relationship with a GO term. So I made sure to filter the "Significant" column in results table so that it was at least equal to 1 or greater. 

To ensure you did it right, calculate the number of significant DGEs for a given time point (i.e. for hour 1 upregulated padj < 0.05 and L2FC>0, it is 26 genes), none of the values in the "significant" column should be greater than 26. Annotated column could be in the thousands, based on the length of genes in the gene universe.
