---
layout: post
title: Fv/Fm GLMMs
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


Moving forward, I'm going to collate notes from [Cunning et al 2021](https://royalsocietypublishing.org/doi/10.1098/rspb.2021.1613), [Klepac and Barshis 2020](https://royalsocietypublishing.org/doi/10.1098/rspb.2020.1379), and [Evensen et al 2021](https://onlinelibrary.wiley.com/doi/full/10.1002/lno.11715) here. 

### Normalization techniques:

1. *Cunning et al. 2021* - **no mention of normalization technique**; there wasn't a component of time in terms of repeated measures of indiviudal corals, so there wasn't a need to normalize to initial Fv/Fm measurement.
2. *Klepac and Barshis 2020* - **Normalized CBASS Fv/Fm values (21–0 h)/0 h** were used for statistical analyses to correct for between ramet variation in starting values. Fv/Fm values measured at the end of each assay were used for plotting for simplicity to allow for easy comparison to previous studies.
3. *Evensen et al. 2021* - **no mention of normalization technique**; each time point was analyzed separately.


### Statistical model used:

1. *Cunning et al. 2021* - Normality of ED50s was tested using Shapiro-Wilk test. Differences in ED50 among nurseries was tested using **Welch ANOVA** after Levene's test indicated unequal variances. Wilcoxon rank-sum tests were used to test for pairwise differences between nurseries, with Bonferroni adjusted p-values and α = 0.01. Relationships between ED50adj and source colony locations, thermal regimes and symbiont abundances were tested using linear models.
2. *Klepac and Barshis 2020* - Differences in normalized Fv/Fm were tested using a **mixed model ANOVA**, where time, a combined origin_destination site variable (owing to the unbalanced design (i.e. not all origins in each destination), and treatment were modelled as fixed factors, and colony identity was nested within experimental tank designation as a random factor. Multiple comparisons across factors and interaction terms were assessed post hoc using general linear hypothesis testing and multiple comparisons (**glht** function; [62]) for linear mixed-effects models, specifying Tukey’s test. To satisfy model assumptions, normality was examined using the shapiro.test and homoscedasticity via the bartlett.test in R, as well as plotting residuals.
3. *Evensen et al. 2021* - Individual response variables were analyzed using linear mixed effects (LME) models, with temperature and experiment as fixed effects, and genet and tank replicate as random effects to account for non-independence of fragments from the same colonies and any potential tank effects. Time points again were analyzed separately. Additionally, LMEs were conducted for each response variable to compare the 27C treatments to the field control and measurements at T0, using the same aforementioned model structure. Models were conducted in the **"lmerTest” package** (Kuznetsova et al. 2017) in R. Model simplification was conducted by backwards elimination of predictor variables (step function), with F-tests used to compare full and reduced models upon removal of the fixed effects. For significant main effects, Tukey’s honestly significant difference (HSD) post hoc pairwise comparisons were conducted using the package **"emmeans”** (Lenth et al. 2020).or all models, distributions of the residuals were plotted to ensure they fit a normal distribution and residuals were plotted against fitted values to confirm that the errors had constant variance. Data were log-transformed to help meet these assumptions when necessary. In order to produce a standardized and comparable proxy to quantitively determine the upper thermal limit of the corals in each experiment, we computed a critical temperature threshold using measurements of Fv/Fm. Measurements of Fv/Fm were fitted to **log-logistic dose–response curves (DRCs) using the package "drc”** (Ritz et al. 2015), with model selection based on Akaike’s Information Criterion (AIC; Table S2) and individual curves fit for each experiment. From these, an "ED50” parameter (effective dose 50) was obtained for each experiment, representing the x-value at the inflection point of the curve (in this case the temperature) where Fv/Fm values in the model fit were 50% lower in comparison to the starting values of the model. This provided a quantitative thermal threshold, designated as the Fv/Fm ED50, for S. pistillata in each experiment.
