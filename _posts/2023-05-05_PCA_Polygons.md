---
layout: post
title: Creating Polygons for PCAs
date: '2023-02-08'
categories: RNA
tags: [RNA]
---

![Screen Shot 2023-05-05 at 11 37 11 AM](https://user-images.githubusercontent.com/56000927/236503510-e0fbf806-a960-4bc4-82f9-513c279b28e3.png)

Figure 1. Principal component analysis of wounding (orange triangle) versus control (blue circle) treatments at hour zero (A), one (B), two (C), and four (D) post-injury. Counts normalized using variance stabilizing transformation (VST) from DESeq2 (Love et al. 2014). 

I want to create polygons with these points so that it's easier to see the clustering. 

Ok that was surprisingly easy, just use geom_polygon() 

![Screen Shot 2023-05-05 at 11 57 14 AM](https://user-images.githubusercontent.com/56000927/236507879-544c798e-8e02-40de-8024-1735236dd8bf.png)

Maybe I can try to add the PERMANOVA results on there too? 

All the code for my PERMANOVAs says that none of the PCAs are significant, but hours 2 and 4 look like they would be significant.

After googling this for a couple hours, I'm going to just leave it as is and not add stats to the PCA because I don't think it's necessary. Also I am not sure if the tests I even did were correct.
