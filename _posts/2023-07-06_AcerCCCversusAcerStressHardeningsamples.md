---
layout: post
title: Acer CCC versus Acer Stress-Hardening Samples
date: '2023-07-06'
categories: coding
tags: [coding, temperaturevariability2023, CCC_ch4]
---

So, the major difference between the CCC samples and the stress-hardening samples is that the first were sequenced using Lexogen 3' QuantSeq, and the second were sequenced following the TagSeq protocol, which is a "specialty" library prep specifically tested and refined with coral samples. So although both projects were sequenced on the same type of machine with the same goal for output (NovaSeq S2 SE100), the results are vastly different. 

I don't think it is based on the extractions because when I nanodropped and Qubitted for both projects, they all had really high yields of RNA. Only a small fraction from each project were checked for RIN scores, so it is possible that the CCC ones had really low RIN scores in comparison. 

The reason I think the CCC ones are worse is when you look at the multiqc of the raw reads, some have less than a million sequences:

<img width="658" alt="Screen Shot 2023-07-06 at 1 49 03 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/b727a8fb-0cbf-424b-bff8-a2facc547e94">

This plot also shows the percent failed versus the sequencing depth. 

<img width="650" alt="Screen Shot 2023-07-06 at 1 48 16 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/735d3883-631c-40e7-a0a1-595b6bc926ef">

Also, when you look at the overrepresented sequences, some are really high which suggests contamination. 

<img width="695" alt="Screen Shot 2023-07-06 at 1 43 57 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/23a3f5a0-a455-43c3-85ac-70adbca14107">

Trimming and polyA tail removal seems to get rid of a lot of base pairs for some samples:

<img width="655" alt="Screen Shot 2023-07-06 at 1 49 43 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/3dbd0e04-c189-43e6-a95a-7da4a4d2c8a6">

Which those happen to be the same samples with high overrepresented sequences and low sequencing depth overall.

Following STAR alignment to the 2019 Acer genome, I seem to get high alignment for some samples but very low alignment for others. 

<img width="672" alt="Screen Shot 2023-07-06 at 1 52 45 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/46d57711-977a-450c-830b-41bd49d76a23">

It looks like a lot of things were unmapped due to being too short, which I think indicates that the sequencing didn't work well or the cDNA library prep didn't work well? 

<img width="686" alt="Screen Shot 2023-07-06 at 1 53 47 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/b4c4150c-d90d-4dac-9592-124356bbd06d">

This ends up being a problem though because when I try to run stringtie and gffcompare (even after filtering out some of the samples -- 1088, 1100, 2264), I get issues with the transcriptome assembly/alignment back to the genome. Specifically, it says that everything aligns perfectly to the original genome, which is not right there is no way. 

<img width="500" alt="Screen Shot 2023-07-06 at 1 55 47 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/969c52d6-0840-4ac2-8cd1-80956d03b9d3">

