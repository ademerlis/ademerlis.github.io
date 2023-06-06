---
layout: post
title: WGCNA for wound healing manuscript
date: '2023-06-06'
categories: RNA, coding, woundhealingPdam
tags: [RNA, coding, woundhealingPdam]
---

I'm currently working through the WGCNA tutorial for the wound healing dataset (following https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html and I specifically did the step-by-step network construction, not the automatic one).

The first thing I noticed in this tutorial is that all the traits that the authors use to correlate with the modules are numeric, but my traits are condition (wounded vs. control) and time of sampling (hour=0,1,2,4). So I feel like these things are abnormal for the input of WGCNA and need to be considered more closely. 

For the purposes of learning how to use WGCNA again, however, I encoded condition as either a 1 or a 2, with 1 = wounded and 2 = control. Hour I kept as 0-4. based on this, I got this dengrodram with modules:

![image](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/d252eba6-6363-4c37-bc97-01f0e9b404f2)

Which I then further merged modules with similar gene expression profiles so I wouldn't have as many. 

Then, the really important figure, is the heatmap of module-trait relationships.
![image](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/c6eac17a-f8b1-4238-9294-ffb1900fad39)

There are more significantly high correlations with hour than there are with condition, but I feel like separating these two traits into columns makes no sense for the research question. It should really be control vs. wounded over time, rather than just time or just injury. 

When I specifically look at the module membership vs. gene significance for the condition column, i do not get any significant correlations.
![image](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/51046980-caff-40f3-a5dc-3fdb841e91f1)

![image](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/4b237474-e69d-4f96-8582-e52c0a7292a5)

![image](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/e54a67d5-1bba-4512-bf3d-cdd2ad1d2eeb)


So I think I need to redo this analysis with the time component in mind.

Also another thing I just noticed -- when going through the SCTLD WGCNA code that Mike Connelly wrote, I realized he had made some more specifications to the network construction part of the code that I had not included in mine. These suggestions are actually also in the WGCNA FAQ section (https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html#:~:text=Data%20analysis%20questions,network%20to%20be%20biologically%20meaningful.)

https://github.com/ademerlis/sctld_transcriptomics_2021/blob/main/ofav_wgcna_updated.Rmd

```{R}
#From Mike Connelly:
# I have chosen the following network construction parameters for the following reasons:
# First, following the recommendations of the WGCNA developers (https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html), a signed network was chosen to be able to detect positive and negative gene correlations, and the biweight midcorrelation was used since it is more robust to outliers. 

#adjacency matrix = how close together the nodes are
adjacency <- adjacency(datExpr_ofav,
      # Network construction arguments:  correlation, adjacency function,  and topological overlap map options
                       corFnc = "bicor", # midweight bi-correlation (more robust to outliers than pearson correlation)
                       power = 20, #based on soft threshold power calculation
                       type = "signed") # signed - it knows down vs. upregulation. more biologically relevant modules. 


#Topological overlap matrix - second step that cleans up the adjacency network. Inferences become more interpretable.
TOM <- TOMsimilarity(adjacency,
                     TOMType = "signed",
                     verbose = 5)
dissTOM <- 1-TOM
# TOM is a measure of similarity. But for clustering you want to see dissimilarity of clusters. So that's why you run dissTOM
rm(adjacency) # may need to delete adjacency, TOM to clear up vector memory
```

So those are things I also need to do when I re-run my analysis (unless the time-series analysis pipeline for WGCNA suggests otherwise)



