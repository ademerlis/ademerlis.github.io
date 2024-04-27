---
layout: post
title: Acer CCC KEGG Pathway analysis issues
date: '2024-04-27'
categories: Analysis, Processing
tags: [Ch4_AcerCCC, coding]
---

I presented my data chapters 2 and 3 to NTK's lab on Friday, and she pointed out something I hadn't thought of - that even though the most significant KEGG pathway in my Acer CCC vs. Nursery analysis was estrogen signaling, the DGEs themselves that had KO terms which annotate to this KEGG pathway could also be involved in many other pathways. Looking at the full list of genes with the KEGG terms [here](https://github.com/ademerlis/AcerCCC/blob/main/gene_expression/results/KEGG_CvN_results_withgenes.csv), you can see that a lot of the genes are involved in many potential pathways. So why did the KEGG analysis give me the estrogen signaling pathway specifically?

When I take all the DGEs (p-adjusted<0.05) from this dataset, and input their KO terms into KEGG mapper, I get this list of KEGG pathways:

<img width="507" alt="Screen Shot 2024-04-27 at 3 45 25 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/94102489-c653-4a45-9735-3932aadc0c2a">

Which is the same as what I get from the results table linked above. Estrogen signaling doesn't have the most KO terms in this list, but it ends up being the most significant. The significance test to obtain p-values is the hypergeometric test, which is essentially testing if the number of genes in a given pathway is more than what would be expected by chance. This number is determined based on the "universe," which is the number of KEGG terms in your dataset from the total number of genes in your dataset. In my dataset, there were 

<img width="559" alt="Screen Shot 2024-04-27 at 3 53 19 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/4010c271-7b81-4cc4-a123-e6d30d5c922b">
