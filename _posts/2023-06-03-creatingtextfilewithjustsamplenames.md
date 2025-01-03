---
layout: post
title: creating .txt file with just sample names using awk/sed
date: '2023-06-03'
categories: [Successful code]
tags: [bash, awk]
---

I want to create a .txt file that just has the sample names listed, no file extensions. This worked for me (ran directly in terminal in folder with fastq sample files):

```{bash}
 awk '{ print FILENAME; nextfile } ' *.gz | cat > samples.txt
 cut -d\. -f1 samples.txt
```
I don't exactly get how it worked but it worked
