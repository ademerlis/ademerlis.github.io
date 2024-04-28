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

Then, when you click on estrogen signaling specifically, these are the genes and associated KEGG terms from my list of DGEs.

<img width="559" alt="Screen Shot 2024-04-27 at 3 53 19 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/4010c271-7b81-4cc4-a123-e6d30d5c922b">

Which is the same as what I get from the results table linked above. Estrogen signaling doesn't have the most KO terms in this list, but it ends up being the most significant. The significance test to obtain p-values is the hypergeometric test, which is essentially testing if the number of genes in a given pathway is more than what would be expected by chance. This number is determined based on the "universe," which is the number of KEGG terms in your dataset from the total number of genes in your dataset. In my dataset, there were 12,344 KO terms from the *A.cervicornis* genome. Then, of my 819 DGEs, 413 of them had KO terms associated with them. The KO terms and individual genes, however, don't exactly point to estrogen signaling. The annotations seem more related to heat-shock and cell signaling (i.e. calmodulin, calcium ion binding). 

So why is it linking these to the estrogen signaling pathway? 

I think perhaps what is happening is that the most significant DGEs in my dataset all happen to be a part of estrogen signaling, but they also are a part of many other pathways. Yes, all of these genes together are significantly enriched in estrogen signaling because they passed the criteria for being more highly expressed than due to chance. But, what gene specifically would point to estrogen signaling. 

I think I need a lot more evidence to support the role of endocrine/hormonal signaling in Coral City Camera corals in comparison to the Nursery corals.

Since I'm running low on time to prepare my defense presentation, I'm going to not go down this path any further at this time. But, I compiled a list of papers to review later and look at which genes are specifically linked to endocrine signaling in these papers. If a lot of my genes are also mentioned, then I think that warrants further discussion about the potential role of this pathway. DO THIS LATER!

papers:
- https://peerj.com/articles/1982/
- https://www.int-res.com/abstracts/meps/v269/p121-129/
- https://www.sciencedirect.com/science/article/abs/pii/S1095643321000167?via%3Dihub
- https://www.sciencedirect.com/science/article/abs/pii/S0048969722071194
- https://www.sciencedirect.com/science/article/abs/pii/S0166445X18305216?casa_token=XQuzXdc9at4AAAAA:QaPEECqFZxjvIwkzwiiaF9pyPDxnmWYYkRDmp2plorfIO4dBzt9s6I7hWEJWlr6h1LLKtqlC
- https://www.int-res.com/abstracts/meps/v269/p121-129/
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8910684/
- https://www.sciencedirect.com/science/article/pii/S0048969721037049#bb0665
- https://ehp.niehs.nih.gov/doi/10.1289/ehp.5233
