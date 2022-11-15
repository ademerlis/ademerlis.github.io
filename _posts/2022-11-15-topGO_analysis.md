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
