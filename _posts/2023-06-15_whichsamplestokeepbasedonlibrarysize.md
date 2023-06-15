---
layout: post
title: Which samples to use based on minimum library size
date: '2023-06-15'
categories: coding
tags: [coding, CCC_ch4, temperaturevariability2023]
---

I have decided which aligner to use, but before I begin that step, I need to determine which samples have made the cut based on the read trimming and QC.

For example, the range of number of reads per sample for the Acer CCC samples varies widely. 

<img width="1103" alt="Screen Shot 2023-06-15 at 8 28 39 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/6a5a71ef-598f-4f7b-9531-6be76f8c90b8">

From memory, I remember something about 5 million total reads being a good baseline. Let me look at some coral papers that have been published recently and see if there is a consensus.

This Illumina article states that the minimum number may be 5 million reads per sample, however it depends on the species too and what has been published for that species. (https://knowledge.illumina.com/library-preparation/rna-library-prep/library-preparation-rna-library-prep-reference_material-list/000001243)

I found this "best practices for RNA-seq" doc from the ENCODE project which I thought was helpful because it gave specific numbers for things: https://www.encodeproject.org/documents/91494746-0ffe-4931-b219-a09802ce1cfa/@@download/attachment/RNA_standards_v1_2011_May.pdf 

Note that it is tailored to the human genome specifically though, so some of the numbers may be more strict than what would be necessary for a coral genome.

**Coral papers**
1. Traylor-Knowles et al. 2021: 2 million reads per sample
2. Young et al. 2020: average of 10 million reads per sample
3. Libro et al. 2013: 4 ± 1 M reads per sample
4. Savary et al 2020: ∼2 billion paired reads (PE-150) with 22,263,882 ± 7,159,826 paired reads/sample (mean ± SD)
5. Helmkampf et al. 2018: 17–42 million


I think I can get away with a 2 million read cutoff. 

For Acer CCC, this removes samples 1088, 1100, and 2264 from the analysis.

For temperature variability, this removes samples Pcli-011, Pcli-045, Pcli-058, Pcli-064, Pcli-066, Pcli-069, Pcli-114, and Pcli-132. 


