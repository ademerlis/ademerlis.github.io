---
layout: post
title: Ch4_AcerCCC updated gff STAR alignment results
date: '2023-07-21'
categories: coding
tags: [coding, Ch4_AcerCCC]
---

The STAR alignment job finished, so I quickly ran multiqc to compare the original alignment, the first round of updated annotations, and the most recent updated gff3 file. What's interesting is all their alignment rates look the same... 

1) **trimmed sequences**
<img width="1077" alt="Screen Shot 2023-07-21 at 11 42 35 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/e3156ada-fe94-4ac1-b529-a68b900b9cb5">

2) **1st round of updated gff3**
<img width="1067" alt="Screen Shot 2023-07-21 at 11 44 31 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/cea8801c-ade6-4861-b8f6-a834d3362378">

3) **most recent gff3 updated file**
<img width="1089" alt="Screen Shot 2023-07-21 at 11 44 51 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/f5d2ebfa-c776-432f-a216-f13182cf5644">

It looks like the only thing that is different is that the multiqc report from the most recent gff3 updated file has "ambiguous features".

<img width="1078" alt="Screen Shot 2023-07-21 at 11 49 02 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/84708126-53dc-455b-b51d-e11b0f388e5a">

So does this mean that the alignment isn't the issue?????? goddamnit. 

I guess what I need to do next is try to extract the gene counts from these and compare those.  
