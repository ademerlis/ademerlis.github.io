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
