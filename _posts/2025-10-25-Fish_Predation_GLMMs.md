---
layout: post
title: Fish Predation GLMMs
date: '2025-10-25'
categories: Analysis, Processing
tags: [fish predation, RTT project, GLMMs]
---

So I've been diving into the RTT project data and trying to find the best pipeline to analyze different response variables. Since there are so many predictor variables, I've been struggling to come up with a way to test all of them and whether they're important for each response variable I'm testing.

Here is the complete list of what I'm working with data-wise:

#### Response variables to test:
- alive vs. dead (binary)
- predation occurrence (binary - fish bites + % live tissue due to predation > 0)
- time to mortality (binary)
- proportion of live tissue (bounded 0-1)
- number of predation events (count)

#### Predictor variables to test:
- Region (categorical)
- Site (categorical)
- Strata (categorical)
- Genotype (categorical)
- Species (categorical)
- Time point (could be continuous or categorical - depending on which column you use)
- Fish abundance (counts)

I also have been struggling with which to use, glm or glmer, and which variables to make fixed effects versus random effects.

When I look at Adam MacAnally's original code, he uses glm exclusively and treats everything as a fixed effect. He uses days as a polynomial (both linear and quadratic) since the data aren't linear. 

But when doing something like repeated measures, you're supposed to include time as a random effect. 

But, do we want to test the significance of time point as a variable and get stats for it? like saying that there is a significant effect of time. If so, then time needs to be a fixed effect. 






### Predation prevalence

##### Model 1 (monitoring interval as random effect)
```{r}
#some variables are nested and need to be controlled for in the model, but can't be fixed effects (model becomes too complicated)

# fixed effects for this test: Species, Region 
# random effects: Source Nursery, Site, STRATA, Monitoring Interval (repeated measures)
# I can't really test genotype because there is very few replicates across sites and regions  

# step-wise add fixed effects and then look at top_models and see which fixed effects are important

# Species and region are important

glmm_global <- glmer(Experienced_Predation_num ~ Species + Region + (1|SourceNursery) + (1|Site) + (1|STRATA) + (1|MI_CORRECTED), family = "binomial", data = outplant_monitoring_days_survival_glm)

summary(glmm_global)
# looking at the variance of random effects, MI_Corrected contributes the overwhelming majority
# source nursery contributes very little to the variance

# fit candidate subsets of predictors
options(na.action = "na.fail")  # required for dredge

model_set <- MuMIn::dredge(glmm_global, rank = "AICc", REML = FALSE)

# top models
top_models <- model_set %>% as_tibble() %>% arrange(AICc) %>% slice(1:8)
# best model with lowest AICc is the model that contains all fixed effects (global model) AICc = 12731.33	

model_no_species <- update(glmm_global, . ~ . - Species)
r2_full <- MuMIn::r.squaredGLMM(glmm_global)[, "R2m"]
r2_no_species <- MuMIn::r.squaredGLMM(model_no_species)[, "R2m"]
species_contrib <- r2_full - r2_no_species
# R2 with both species and region is 0.233
# R2 with just region is 0.223
# Species thus contributes only 0.01 to the explained variance


# residual diagnostics (DHARMA)
sim_res <- simulateResiduals(glmm_global, n = 1000)
plot(sim_res)
testUniformity(sim_res) # fits, not significant
testDispersion(sim_res) # fits, not significant
testZeroInflation(sim_res) # fits, not significant

# make summary table 

# fixed effects
fixed_tbl <- tidy(glmm_global, effects = "fixed", conf.int = TRUE)

# random effects
random_tbl <- tidy(glmm_global, effects = "ran_pars", conf.int = TRUE)

# combine
full_tbl <- bind_rows(
  fixed_tbl %>% mutate(effect_type = "fixed"),
  random_tbl %>% mutate(effect_type = "random")
)

# Perform pairwise comparisons to look at the impact of treatments within MCAV
pairwise_emmeans <- emmeans(glmm_global, pairwise ~ Species*Region, adjust = "tukey")

Anova(glmm_global, type = "III")
# significant effect of Species and Region

# View the pairwise comparisons results
print(pairwise_emmeans$contrasts)
```


