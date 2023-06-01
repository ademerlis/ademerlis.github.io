---
layout: post
title: initial FastQC results for Acer CCC samples (untrimmed)
date: '2023-06-01'
categories: coding
tags: [coding]
---

So I got the multiqc report for the CCC Acer samples to work, and it seems like they are all over the place.

**Note: these files have not been trimmed in any way.**

I need to double-check they are 3' RNA-Seq and what kind of sequencing they were done on. (going to look at Natalia's presentation)

It is 3’ Quantseq RNA sequencing

Also a note from Natalia: the Novaseq S2 brings a lot of "artifacts" that need to be trimmed out. PolyA tail can be really long, even post-trimming, so that compromises alignment.

She recommends using TrimGalore and use a code specifically to remove artifacts from NovaSeq + PolyA tails.

So, let's try trimming them and then see if the FastQC results below become better (before we try to interpret these results below).

<img width="1061" alt="Screen Shot 2023-06-01 at 12 45 28 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/dd49481d-8f43-464a-b684-e3eeb46adcbc">

<img width="1086" alt="Screen Shot 2023-06-01 at 12 45 50 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/357aff7b-740f-4880-9f1c-253e102a64af">

<img width="1084" alt="Screen Shot 2023-06-01 at 12 46 01 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/99024ac9-6881-4b0a-b76f-314389db80a4">

<img width="1089" alt="Screen Shot 2023-06-01 at 12 46 15 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/f57b809b-79f0-4534-ab29-6f3fce34a615">

<img width="1081" alt="Screen Shot 2023-06-01 at 12 46 27 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/ea2d388b-db01-4873-8806-95defb212ded">

<img width="1082" alt="Screen Shot 2023-06-01 at 12 46 38 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/8b3de1aa-663d-49b8-bb89-05d81f3086c6">

<img width="1083" alt="Screen Shot 2023-06-01 at 12 46 54 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/b749d35e-a674-4001-94ba-7a551699461c">

<img width="1078" alt="Screen Shot 2023-06-01 at 12 47 13 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/10d5d056-c6ab-4fb9-82b3-558f22c28e59">

<img width="1068" alt="Screen Shot 2023-06-01 at 12 47 26 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/dee4c22a-c9b0-48a3-a307-ee0c67c8aed0">



