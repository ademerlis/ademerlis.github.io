---
layout: post
title: some wins against Pegasus
date: '2023-06-07'
categories: [Successful code]
tags: [Pegasus, bash, stress-hardening, AcerCCC]
---

After meeting with Anthony Bonacolta, I finally got some help as to why my job submissions to Pegasus weren't working. First, the alias "compute" that I created includes "bsub -P and_transcriptomics" and I was trying to do compute .sh every time I submitted a job. When I did it this way, it also wasn't recognizing when i put something in the bigmem queue and would put it in the general queue. 

Anthony said that because in my .sh file, I had already specified the -P flag in the # at the top, I don't need to do it again when I submit the job.

Instead, just run:
```{bash}
bsub < job.sh
```
And this has worked for me every time so far! 

So now for the CCC samples, I've trimmed the adapters and short reads, and then I wrote a second job to trim the polyA tails (because in TrimGalore you need to trim things sequentially, idk why). I've also run fastqc and multiqc on the trimmed reads, and I will run that again on the trimmed_trimmed reads (w/o polyA tails hopefully).

For the stress-hardening samples, I ran fastqc and multiqc on the raw reads, and I am currently running TrimGalore to remove the adapters and short reads. Then I'll do the same thing with the polyA tail code and then run fastqc and multiqc again. 
