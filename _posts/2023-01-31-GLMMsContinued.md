---
layout: post
title: More notes on GLMMs (with Dr. Kevin Wong!)
date: '2023-01-31'
categories: [Trialing code, Statistics]
tags: [GLMM, Linear Mixed Models, PAM, Stress-Hardening, glmmTMB, Kevin]
---

We looked through my PAM_stresshardening.Rmd code that I wrote for running glmmTMB (fvfm ~ Treatment + (1|Colony) + (1|Tank/Date), family = list(family="beta", link="logit")).

After running it, there seems to be a couple red flags. One is that when we run the summary of the model, the fixed effect of Treatment is written out as "TreatmentVariable" in the results table. That could just be how the model output formats it, but it's kind of strange. So you need to look into that further to make sure it is using both variable and control data.

Then, when looking at the formula I wrote, Kevin actually thought that Date shouldn't be nested within Tank, because part of the research question is how Fv/Fm changed over time. So the formula should actually be fvfm ~ Treatment * Date + (1|Colony) + (1|Tank), so Colony/Genotype is the random effect and Date is included as a fixed effect. 

The next thing that was missing was that I only ran one model, the full model with all the explanatory terms. I need to run each iteration of the model with different levels to see which model best fits the data. To do this, I run each iteration:
1. Treatment*Date + (1|Colony) + (1|Tank)
2. Treatment*Date + (1|Colony)
3. Treatment*Date + (1|Tank)
Then, compare the AICs and the likelihood ratio test p-values to see which model is best. If the complex model is not significantly different from the less complex model, then go with the less complex model.

The last red flag that I need to look into with my model is to check assumptions. For proportion data models, the assumption is that there is overdispersion (variance > mean). I need to check that.

One other thing I could do is normalize my Fv/Fm data so that it starts at 1 and then decreases based on the previous metric (so normalize it to the initial value). That could potentially make the data fit better (when running qqnorm on the residuals).
