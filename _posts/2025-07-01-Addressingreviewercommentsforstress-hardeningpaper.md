---
layout: post
title: Addressing reviewer comments for stress-hardening paper
date: '2025-07-01'
categories: Analysis
tags: [Coding, Ch2_tempvariability, DESeq2, PCoA]
---

A reviewer had a concern about the intra-genotypic variation in the Acer and Pcli stress-hardening manuscript that I submitted to Ecology and Evolution (see figure below).

![image](https://github.com/user-attachments/assets/3dc8fb77-0745-4bcc-b599-846f1012aa7a)

At first I thought they were referring to genotypic variation in general, which I did report p-values for because genotype was a factor in the PERMANOVA. But, I think their concern is specifically for within a treatment, the samples of a given genotype are far apart on the PCoA.

So then I asked ChatGPT if there was a test to do that, and of course it had an answer - use "PERMDISP" to assess "multivariate homogeneity of group dispersions." In the vegan package, there is a function called "betadisper" and you can select the grouping variable you want to assess. It's kind of like the Levene's test for multivariate analyses.

So I ran this for genotype for Acer host gene expression first:

First run the PCoA function to make sure it's making the same PCoA as the manuscript and also to make sure the variables are all loaded.

```{r}
load("Rdata_files/host/initial_fullddsdesigncountsVsdcounts.RData")
ad.pcoa=pcoa(dist(t(assay(Vsd)),method="manhattan")/1000)
scores=ad.pcoa$vectors
ad.pcoa$values # the % variation explained by each axis is the Relative_eig column
conditions=design

scale_color_manual(values = c("grey", "#FF3333", "#00CCCC"))
#code by group
pdf("PCoA_Treatment_Genotype.pdf")
plot(scores[,1], scores[,2], pch=c(15,17,25)[as.numeric(as.factor(conditions$Genotype))], col=c("#00CCCC", "grey","#FF3333")[as.numeric(as.factor(conditions$Treatment))], xlab="PC1: 20%", ylab="PC2: 11%", main="Treatment")
ordispider(scores[,1:2], conditions$Treatment,col=c("#00CCCC", "grey", "#FF3333"))
legend("topright", legend=c("Control", "Initial", "Variable"), fill = c("#00CCCC", "grey", "#FF3333"), bty="n")
legend("topleft", legend=c("BC-8b", "MB-B", "SI-C"), pch=c(15,17,25), bty="n")

#analysis of variance 
ad=adonis2(t(assay(Vsd))~Genotype+Treatment,data=conditions,method="manhattan", permutations=1e6)
ad
as.data.frame(ad) %>% 
  write_csv("PERMANOVA_Acer.csv")

#analysis of dispersion (within-group variation)

# Test for differences in dispersion by genotype
dist_mat <- (dist(t(assay(Vsd)),method="manhattan")/1000)
dispersion_genotype <- betadisper(dist_mat, group = conditions$Genotype)
anova(dispersion_genotype)
permutest(dispersion_genotype)

# You can also visualize
plot(dispersion_genotype)

TukeyHSD(dispersion_genotype)

# Test for differences in dispersion by genotype within a treatment
conditions$Geno_Treatment <- paste(conditions$Genotype, conditions$Treatment, sep="_")
dispersion_geno_treatment <- betadisper(dist_mat, group = conditions$Geno_Treatment)
anova(dispersion_geno_treatment)
permutest(dispersion_geno_treatment)

# You can also visualize
plot(dispersion_geno_treatment)

TukeyHSD(dispersion_geno_treatment)

```

When I just run genotype and don't separate out treatments, it is not significant.

![Screenshot 2025-07-01 at 9 01 33 AM](https://github.com/user-attachments/assets/247e38f1-fa60-44b3-a257-219945f23ad6)

![Screenshot 2025-07-01 at 9 02 01 AM](https://github.com/user-attachments/assets/5defda68-220f-4fb6-b236-80c4da2bc611)

But then when I run the Tukey, while it is still not significant, you can see that is comparing the different groups - the different genotypes. So it isn't actually assessing within-genotype differences, it's comparing the dispersion between genotypes.

![Screenshot 2025-07-01 at 9 02 14 AM](https://github.com/user-attachments/assets/8a3e08d4-2471-4d00-a232-e59b8f3045fe)

That's not what I want. And I guess that's also the point of the Levene's test, is that you're comparing the homogeneity of variation between groups - is the dispersion of one sample group similar to the dispersion of the other? 

What the reviewer is asking is, within a sample group, why are the samples so different? Is that statistically significant?

So then I ask ChatGPT again, and this time it tells me to use "anosim" for within-group structure, or substructure, within a grouping variable (in this case, treatment). 

But I don't think that's what I want either? Because that sounds like then it is testing whether treatment is significant within a given genotype.

The last option that ChatGPT gives me is to do a custom permutation test to evaluate whether a genotype is particularly dispersed, compared to a null expectation. I am wary of trying that because the code it provided me creates a random gene set based on the genes in my data set, and then compares dispersion of my genotypes to the random set. I could do that I guess, but 1) it's time consuming and 2) I don't think that will actually prove anything.

Essentially, there is going to be inherent variation between samples. I don't know how to statistially show that, but that's the point of having biological replicates so that you get a mean across all the samples. I think if someone was to even sample the same coral at the same time point, you would still get some sort of difference between samples on a PCoA. Obviously you would hope that it wouldn't be more different than the treatment effect, but there will still be some difference. The samples won't be right on top of each other.

I think I will attach the PCA graph as another visualization, which shows much less spread and may be what the reviewer is more used to seeing. 


![Screenshot 2025-07-01 at 9 13 04 AM](https://github.com/user-attachments/assets/3195827e-e1a4-467a-936a-f8a852befcee)
