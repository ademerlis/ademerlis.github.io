---
layout: post
title: WGCNA for wound healing manuscript
date: '2023-06-06'
categories: [WGCNA]
tags: [WGCNA, wound healing]
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

Direct quotes from FAQ that are relevant:
"The default correlation method in all functions in WGCNA is standard Pearson correlation. In general, unless there is good reason to believe that there are no outlier measurements, we recommend (and use ourselves) the **biweight mid-correlation** as a robust alternative. This is implemented in WGCNA function bicor. Many WGCNA functions take the argument corFnc that allows one to specify an alternative correlation function to the standard cor and bicor is one option. Additional arguments to the correlation function can be specified using the argument corOptions (depending on function, this argument may require one of two alternate forms, please see the help for each function for details). In certain functions, notably the of the blockwise family, correlation function cannot be specified directly as a function; rather, one must use the argument corType to specify either Pearson or biweight mid-correlation.

Some cautionary notes regarding the use of bicor:

1. **Restricting the number of excluded outliers:** argument maxPOutliers. The default version of the biweight mid-correlation, described in Langfelder and Horvath (2011) (link to article), can produce unwanted results when the data have a bi-modal distribution (e.g., when a gene expression depends heavily on a binary variable such as disease status or genotype) or when one of the variables entering the correlation is itself binary (or ordinal). For this reason, we strongly recommend using the argument **maxPOutliers = 0.05 or 0.10** whenever the biweight midcorrelation is used. This argument essentially forces bicor to never regard more than the specified proportion of samples as outliers.
2. **Dealing with binary data.** When relating high-throughput data x to binary variable y such as sample traits, one can use argument robustY = FALSE to turn off the robust treatment for the y argment of bicor. This results in a hybrid robust-Pearson correlation as described in Langfelder and Horvath (2011). The hybrid correlation can also be used when one of the inputs is numeric but known to not have any outliers."

I added the maxPOutliers argument to the adjacency calculation, but I'm sure if I should also add the robustY = FALSE. Technically the condition trait would be a binary trait (0 = control, 1 = wounded), so I think I should add it.

Then again these may all change once I look into time-series with WGCNA.

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


Ok I found another tutorial that has helped me figure out how to create the traits matrix for correlating with the modules. This one comes from UT Austin and it's a coral larval development dataset (https://wikis.utexas.edu/display/bioiteam/Clustering+using+WGCNA)

Their traits matrix that they use to correlate has binary (0-1) encoding for their 4 replicates per time point. It looks like this:

<img width="844" alt="Screen Shot 2023-06-06 at 10 18 28 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/6495eb19-425d-44a7-a998-5031dc12d212">

So I need to make a table that has the same thing, with wounded replicates and control replicates as the rows, and then the time points (hour 0, 1, 2, 4) as the columns.

I initially made the metadata file with just the time points as the columns, but then when I created the heatmap I realized that it was combined control and wounded coral samples within that. So, I needed to create sort of an interaction term for each column, i.e. hour0_control vs. hour0_wounded. Then I encoded the binary 0 or 1 based on whether the sample fit that criteria or not.

Ok, after following all the above changes and making a new metadata file, I got this:

<img width="1197" alt="Screen Shot 2023-06-06 at 10 59 58 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/a8407339-e316-4f06-98a8-52080e31c2d7">

What is interesting is that modules seem to correlate more based on time rather than treatment, but that there are some that are opposite between treatments.

