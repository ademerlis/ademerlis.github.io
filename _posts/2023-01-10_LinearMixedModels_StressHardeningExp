---
layout: post
title: Figuring out Linear Mixed Models using Stress-Hardening Data
date: '2023-01-10'
categories: Coding
tags: [Coding]
---

I am revisiting the 2022 stress-hardneing experiment to finally tackle the unknown, which is how to incorporate so many different variables into one statistical test (tank replicates, biological replicates, species, genotypes, treatments, time points).

Dr. Kevin Wong generously sent me a PPT he made to go over these concepts.

First, the rules of linear models. 
Assumptions that must be met:
1. Linearity of relationship between predictors and responses.
2. Normality of residuals.
3. Homogeneity of variance.
4. Independence of residual errors.

Mixed effect models: general or generalized linear models that allow for both fixed and random effects. 
GLMMs are useful for:
1. data with sources of random variability (i.e. tank, genotype)
2. Data with non-independence (tank effects, repeated measures, multilevel or hierarchical designs)
3. Data with non-normal distributions

A fixed effect is a parameter that does NOT vary across observations (i.e. experimental treatment). Time is a fixed effect in a time-series analysis. 

A random effect is a source of variability unexplained by a fixed effect, or data collected across incomplete levels of effect. This would include things like repeated measures (individual sampled multiple times, tank effect), or site, depth, or genotype.

The random intercept is an important part of the linear model, because it allows the model to estimate a diference response for each level of random intercept. (example: growth may be higher in some tanks due to light and flow, so add (1|tank) to your GLMM formula)

Nested random intercepts are important when a lower level factor appears only within a particular level of upper level factor. Example: a coral stays in the same tank for the whole experiment, and is sampled multiple times = (1|tank/sample)

A crossed random intercept is when a lower level factor appears within a multiple level of upper level factor. Example: a coral genotype that is sampled multiple times is distributed across all tanks (1|tank/sample) + (1|genotype).

Random slope allows the model to estimate different effect of fixed effects for each level of random slope. For example, the growth of the coral of small versus large corals is differentially dependent on feeding treatment = (1+size)

When running a model, we need to pick the right distribution. Otherwise, our conclusions will be flawed if the data does not meet the assumptions of the model.

Types of data for use in GLMMs:
1. Measurement (continuous data) --> normal/Gaussian distribution (example given is coral tissue properties between treatments)
2. Binomial/Logistic (0,1; success or failure) --> binomial distribution (example given is survival based on a treatment over time)
3. Count (whole integers, non-continuous) --> Poisson distribution (example given is number of polyps that expand per treatment)
4. Percentage/Proportion --> beta or binomial distribution; the family="beta", link="logit" function in code. (example given is probability density of survival in different treatments)
5. Time to failure --> Kaplan-Meier, Cox Proportional Hazards (example is probability of survivorshop over time in different treatments)
