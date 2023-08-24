---
layout: post
title: DESeq2 design formula
date: '2023-08-23'
categories: coding
tags: [coding]
---

I have started working through [Dr. Michael Studivan's DESeq code](https://github.com/mstudiva/SCTLD-intervention-transcriptomics/blob/main/code/intervention/deseq2_intervention_host.R), but I am having trouble figuring out how best to include all of my explanatory variables in the DESeq model. In his study, he has time and treatment as explanatory variables. He created a column that combined the two: time.treatment. See screenshot of data table with design:

<img width="594" alt="Screen Shot 2023-08-24 at 9 13 31 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/9e167d43-4bce-4e77-894c-982139b8e9e4">

At first, I was confused why he didn't use genotype as an explanatory variable as well, but then I realized that he actually had two different experimental designs going on in the same paper. So this one linked above is called "deseq2_intervention_host", which is the intervention experiment. In this one, they sampled colonies at two different time points in the field. Because it can be assumed that each colony a certain distance apart is a different genotype, essentially all the colonies sampled are different genotypes and there are no replicates within each genotype. (I'm not sure if they did genotype typing). The experimental design I **Should** be looking at for my own analysis is the "transmission_host" files, because that was his ex-situ ERL experiment where he had genotypic replication. Here is the github link: https://github.com/mstudiva/SCTLD-intervention-transcriptomics/blob/main/code/transmission/deseq2_transmission_mcav_host.R 

When I look at the design formula in this one, he does include in the dds model: ~genotype+fate

The whole reason I'm looking into this is because so far, I've been running my DESeq2 model as: group (time_point.treatment) + Genotype. But I'm not sure if this is the correct way to input the explanatory variables. Or the correct order for them. 

I found this 2020 publication ([Law et al. 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/)) that specifically is a guide for how best to design matrices for gene expression experiments. 

Some definitions:
- covariates versus factors: covariates are "numerical values that are quantitative measurements associated with samples in the experiment." Factors are "categorical variables or classifiers associated with samples in the experiment. They are often separated into those that are of a biological nature (e.g. disease status, genotype, treatment, cell-type) and those that are of a technical nature (e.g. experiment time, sample batch, handling technician, sequencing lane)."
- So in my experiment, I don't think I have any covariates for gene expression specifically. All my variables are factors (time, treatment, genotype)

From reading the above paper, because Genotype is not specifically my term of interest, but rather the treatment over time is my specific interest, I would NOT add genotype as part of the grouped variable. The design I have used group(time_point.treatment) + Genotype is correct. Yay!
