---
layout: post
title: Next steps with PSTR RT transcriptomics
date: '2026-01-08'
categories: Analysis
tags: [reciprocal transplant, coding, Pegasus, rRNA, de novo transcriptome assembly]
---

Some of the trimmed sequence files (<2GB) were run on sortmeRNA using --blast, and the others were run more recently (using [the most up to date sortmeRNA code](https://github.com/ademerlis/reciprocaltransplant/blob/main/gene_expression/6_sortmerna_redo.sh)), which did not include --blast.

This resulted in differing files (no singletons) -- asking chatGPT, they said that sortmeRNA doesn't change the filtering parameters for rRNA vs non_rRNA, it just doesn't pull out paired vs singletons.

So in essence, I can use the "old" files and the "new" files together to assemble a transcriptome using trinity.

creating Trinity environment in conda: 
`conda create -n trinity -c conda-forge -c bioconda trinity`

moved all the _non_rRNA and _rRNA files up, and then identified the _out_paired from the first round and renamed those to non_rRNA and moved those too

now identify samples to use for transcriptome assembly (i selected a subset across genotypes and time points, and then looked for ones that had fwd and rev sortmerna files > 1 GB but less than 2 GB).

look for samples from sortmerna with < 1 GB -- rerun those (aka the one at 17 k)

samples to rerun (all at 17kb):
- Pstr_Jul2023_006
- Pstr_Nov2023_180
- Pstr_Nov2023_181
- Pstr_Nov2023_190
- Pstr_Nov2023_191
- Pstr_Nov2023_195
- Pstr_Nov2023_196
