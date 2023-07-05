---
layout: post
title: Troubleshooting STAR alignment
date: '2023-07-03'
categories: coding
tags: [coding, CCC_ch4]
---

This morning to try to troubleshoot the "missing gene IDs" issue with my STAR alignment, I followed one suggestion of ChatGPT to use the [Integrative Genomics Viewer app](https://software.broadinstitute.org/software/igv/) to visualize the alignment of one of the sample BAM files with the Acer genome. 

To do this, I downloaded the IGV GUI to my local computer, then also installed samtools into my bash_profile because I needed to create an index for the .bam file to load it in IGV. To run IGV, you input the .fasta genome file, then your sample .bam file and the index file. For the sample I used, it looked like this:

![igv_snapshot](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/827788f0-c898-44fc-9834-4ec7dffd1cee)

It doesn't look like any of the IGV examples, so that makes me think something is wrong.

I see two options for moving forward (while waiting to hear back about help from a friend): 
1. try pipeline with other Acer samples and see if those are better-quality and thus alignment would work better
2. try something like Galaxy that does most of the coding work for you

I tried using Galaxy, because they have STAR for RNA-seq on there, and the upload rate for just one trimmed.fastq.gz file was so slow that I gave up. Today (July 5) I'm going to try to run FeatureCounts on the aligned reads and see if it works anyways. I'm also going to try to run stringtie with the stress-hardening exp Acer samples and see if those work better.


