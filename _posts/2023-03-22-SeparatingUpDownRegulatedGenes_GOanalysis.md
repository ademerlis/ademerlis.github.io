---
layout: post
title: Separating Up- and Down-Regulated Genes for GO Analysis?
date: '2023-02-08'
categories: RNA Sequencing
tags: [RNA, Sequencing]
---

I have returned to the GO analysis for each hour, as I was talking with Kevin and he mentioned to separate the up and down-regulated gene lists so that you know which pathways are being upregulated versus downregulated. But then I ended up in a deep-dive of GO enrichment analysis, and whether it's appropriate to split up the DEG list by directionality. Turns out this is a contested debate and ultimately depends on your question (I think). 

The main issue is that any enrichment analysis is essentially comparing the list of DEGs to the gene universe you provide, and seeing whether there is more enrichment in the DEG list. That is inherently biased for several reasons... (read https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0761-7) 

But, putting that to the side, there is also the question of whether or not a down-regulated gene is actually contributing to overall higher enrichment of a pathway (maybe the down-regulated gene is an inhibitor, for example). 

See these articles for more discussions/debates:
https://www.researchgate.net/post/While-doing-Gene-Set-Enrichment-Analysis-should-one-separate-up-and-down-regulated-genes
https://www.biostars.org/p/202681/
https://www.biostars.org/p/100123/
https://royalsocietypublishing.org/doi/full/10.1098/rsif.2013.0950?rss=1
https://www.biostars.org/p/155501/
https://peerj.com/articles/159/

Based on these articles, I am going to try to separate them and see if the pathways become clearer, because right now when all DEGs are grouped together, it is difficult to interpret biological meaning. 

I also saw a suggestion about doing WGCNA first and then doing GO analysis on the modules of highly correlated genes. 