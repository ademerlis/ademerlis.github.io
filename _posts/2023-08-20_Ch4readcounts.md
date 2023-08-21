---
layout: post
title: Ch4_AcerCCC read counts
date: '2023-08-20'
categories: coding
tags: [coding, Ch4_AcerCCC]
---

I'm writing up a summary for the Ch4 samples now, and I want to generate a summary table with the average and standard deviation of raw reads, trimmed reads, and percent alignment rates.

```{r}
library(tidyverse)

sample_metadata <- read_csv("../input_files/sample_metadata.csv")

sample_metadata %>% 
  select(Location, Genotype, `Tube No.`, `M Seqs (Raw Reads)`, `M Seqs (post-polyA trimming)`, `% Aligned`, `M Aligned`) %>% 
  mutate(percent_alignment = as.numeric(gsub("%","",`% Aligned`))) %>% 
  group_by(Location) %>% 
  drop_na() %>% 
  summarise(count = n(),
            total_raw_reads = sum(`M Seqs (Raw Reads)`), 
            total_trimmed_reads = sum(`M Seqs (post-polyA trimming)`), 
            average_rawreads_persample = mean(`M Seqs (Raw Reads)`), stdev_rawreads_persample = sd(`M Seqs (Raw Reads)`), 
            average_trimmedreads_persample = mean(`M Seqs (post-polyA trimming)`),
            stdev_trimmedreads_persample = sd(`M Seqs (post-polyA trimming)`),
            average_alignment_rate = mean(percent_alignment),
            stdev_alignmentrate = sd(percent_alignment))  
```

Of 17 samples were sequenced, three were removed from downstream analyses due to poor quality and low alignment rates. Of the CCC samples (N=9), there were 164.4 million raw reads and 140.8 million reads after trimming. The average (± standard deviation) raw reads per sample was 18.27 (± 5.34) million reads. The average (± standard deviation) trimmed reads per sample was 15.64 (± 5.55) million reads. Average (± standard deviation) alignment rate was 40.79% (± 22.2%). Of the nursery samples (N=7), there were 126.4 million raw reads and 112.1 million reads after trimming. The average (± standard deviation) raw reads per sample was 18.06 (± 6.29) million reads. The average (± standard deviation) trimmed reads per sample was 16.01 (± 7.08) million reads. Average (± standard deviation) alignment rate was 47.59% (± 20.34%).
