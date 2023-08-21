---
layout: post
title: meeting with Michael
date: '2023-08-21'
categories: coding
tags: [coding]
---

I met with Dr. Michael Studivan last week to discuss plans for DNA and RNA extractions of the Ch 3 Pstr samples. During this meeting, we talked about the results from the Ch 2 RNA-seq data, and he was concerned with how high the number of passed reads there were post-filtering step (see [the multiqc report](https://github.com/ademerlis/temperaturevariability2023/tree/main/gene_expression/bioinformatics/QC#2-multiqc-reports-of-trimmed-reads)). He said that when he has run TagSeq analysis, he usually get around 60% loss of reads when going from raw reads to trimmed reads. Whereas for mine, I got >95% retention of reads which passed the fastp filter. So what's happening here? Did fastp not work? 

Michael said that TagSeq is know to have a lot of PCR duplicates, and that this is important to filter out so there is only one copy per transcript.

Maybe that is what is happening here in the "sequence duplication levels" figure from multiqc? 

This first graph is from the raw reads:

<img width="816" alt="Screen Shot 2023-08-21 at 10 44 36 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/ec1a0d48-633e-4e48-aa90-867d70e441ba">

There is no graph like that for the fastp-trimmed reads. But, in the first table of general statistics, there is a column for percent duplication. I sorted it by highest to lowest, and it looks like all the Pcli samples have > 60% duplication. But, that doesn't seem to have been trimmed, because the %PF (passed filter) is ~97%. 

<img width="826" alt="Screen Shot 2023-08-21 at 10 47 14 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/3909acc9-bbde-47cf-ac5f-814441476857">

For Acer, the duplication rates are "lower" (30-50%), but the %PF is still really high.

<img width="810" alt="Screen Shot 2023-08-21 at 10 48 24 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/1bd918db-093c-4fd9-b035-7a9fa53c9c9e">

Going back through my [Ch4_AcerCCC multiqc reports](), it looks like TrimGalore/cutadapt 
