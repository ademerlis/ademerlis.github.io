---
layout: post
title: STAR align updated to counts matrix
date: '2023-07-21'
categories: coding
tags: [coding, CCC_ch4]
---

So following the updated STAR alignment with the updated gff3 file for the Acer CCC samples, I re-tried the code from [Natalia](https://github.com/ademerlis/ademerlis.github.io/blob/master/_posts/2023-06-29_STARoutputtoreadcounts.md) and got different results which I think means it worked!

Before:

<img width="1070" alt="Screen Shot 2023-07-21 at 6 06 27 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/2d8251aa-3156-42ce-b669-658608ca5fb6">

After:

<img width="296" alt="Screen Shot 2023-07-21 at 6 06 48 PM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/709233f1-62ba-4b6d-9973-e92c77417b63">

So that means the gene IDs worked, which is good. I haven't tried the stringtie pipeline to see if the MSTRG things come up. I could just continue on with Natalia's pipeline, at least for the Acer ones and see if I get DEGs with that.

I see now though that I also did the stringtie for Acer CCC on the original samples (not with this newest updated gff3), so maybe I could try that and then compare the gene matrices. 
