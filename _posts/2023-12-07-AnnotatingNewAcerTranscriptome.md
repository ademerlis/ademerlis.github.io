---
layout: post
title: Annotating new Acer transcriptome
date: '2023-12-07'
categories: Coding
tags:
  - Coding
published: true
---

I'm going to be following [Michael Studivan's code](https://github.com/mstudiva/Acropora-cervicornis-annotated-transcriptome/blob/main/tagSeq_TranscriptomeAnnotation_README.txt) for annotating the new Acer transcriptome. 

Step 1 is installing bioperl. I had trouble doing this, but I didn't realize I could activate it as a conda environment. It looks like it's working so far.

```{bash}
# To install Bioperl as a conda environment
conda create -y -n bioperl perl-bioperl

conda activate bioperl


```
