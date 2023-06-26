---
layout: post
title: Aligning raw reads versus trimmed reads results for Acer CCC samples
date: '2023-06-26'
categories: coding
tags: [coding, CCC_ch4]
---

The results of the STAR alignment for the trimmed reads were a little strange looking (especially in comparison to the stress-hardening samples, those look way better), so I tried aligning the raw reads and seeing if I got better results (higher alignment and less unmapped reads). Interestingly, it looks like trimming improved alignment rates (trimmed samples on left, raw reads on right):

<img width="1412" alt="Screen Shot 2023-06-26 at 9 20 09 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/d0118a1a-aeb4-4cc4-aa8c-eb835652c4c8">

In comparison, here is what the stress-hardening alignment looks like (for the Acer samples specifically):

<img width="1097" alt="Screen Shot 2023-06-26 at 9 35 52 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/eeb0cacd-281a-4896-9845-3f95e00cf64d">

In terms of what I mean by "better results", I mean higher percentage alignment and overall higher read alignment in terms of number of bases. I'm not sure what a good number is for these but I think it is something like >50% alignment. Idk though.

Also, I'm not sure what a good alignment score is. For the stress-hardening samples, this is what the multiqc report looks like:

<img width="1097" alt="Screen Shot 2023-06-26 at 9 37 55 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/2ce15468-ce72-41e0-8ed1-562030fdb6f1">

versus the Acer CCC samples:

<img width="1129" alt="Screen Shot 2023-06-26 at 9 48 47 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/4df39e5c-88e2-4818-983b-3b9ad75c0a20">

I think overall I can move forward with samples that have greater than 30% alignment for the Acer CCC ones, and move on with all the stress-hardening samples and see how the gene counts come out. Then I need to try to align the unmapped reads to the symbiont genome (Symbiodinium fitti) to see if those get high alignment. 

Also, I thought I wrote this somewhere already but I can't find it, but it is important to note that the Acer CCC and the Acer stress-hardening were sequenced differently.

Acer CCC = 3’ Quantseq NovaSeq S2

Acer stress-hardening = ... 

