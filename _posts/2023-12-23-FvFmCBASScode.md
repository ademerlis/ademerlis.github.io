---
layout: post
title: Fv/Fm CBASS code
date: '2023-12-23'
categories: Coding, IPAM
tags: [Coding, IPAM]
---

I am revisiting my Fv/Fm code for my Chapter 2 stress-hardening analysis. When I first started this, I tried to follow code from [Cunning et al. 2021](https://github.com/jrcunning/CBASS_FL_Acer), who assessed A.cervicornis outplants using CBASS across Florida. However, I found other publications that applied different ways of analyzing this data that were more straightforward (I think because Cunning et al. 2021 had a lot more variation to account for due to the different sites across space and time, and also the authors added corrections for light parameters in each CBASS tank and for use of a Diving PAM). I started playing around with other codes and eventually got really deep in the weeds of this.

My overall goals with my Fv/Fm data are to:
1. evaulate whether the variable temperature treatment impacted Fv/Fm in Acer and Pcli
2. evaluate thermal tolerance of temperature-treated versus untreated corals using CBASS and dose-response curve (DRC) modeliong

With that being stated, my move forward should be simple. 

I think a couple obstacles that held me up were:
1. Do I / how do I normalize Fv/Fm data? Based on initial measurement (pre-treatment)? Based on pre-CBASS measurement? Do I incorporate number of days and create a rate of change and then statistically compare that?
2. Which statistical model do I use? I want to incorporate fixed (treatment, CBASS temperature) and random (genotype, tank) effects, so Generalized Linear Mixed Models is appropriate. However, I was stuck on using glmmTMB with the "beta" family distribution and the "logit" function. The major assumption of this type of distribution is that data are bound by 0 < x < 1. I originally thought fv/fm data would fit this best because it is a proportion that the IPAM is measuring. But, some of the values equal zero, and when I start doing time normalizations I get values greater than 1 or less than zero. So, for this reason I no longer think glmmTMB is appropriate. The next best GLMM would be the classic continuous response model using lme4::lmer. This is consistent with code I've seen from Klepac and Evensen.

See my previous GitHub posts for more confusions:
1. https://github.com/ademerlis/ademerlis.github.io/blob/main/_posts/2023-11-06-TroubleWithFvFmData.md
2. https://github.com/ademerlis/ademerlis.github.io/blob/main/_posts/2023-01-31-GLMMsContinued.md
3. https://github.com/ademerlis/ademerlis.github.io/blob/main/_posts/2023-01-29-GLMMs.md
4. https://github.com/ademerlis/ademerlis.github.io/blob/main/_posts/2023-01-10-LinearMixedModels_StressHardeningExp.md
5. 
