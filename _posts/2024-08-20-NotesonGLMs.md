---
layout: post
title: Notes on GLMs
date: '2024-08-20'
categories: Stats
tags: [stats]
---

After reviewing some great tutorials, I think I understand the best steps for generalized linear models with fixed and random effects (I think).

Tutorials I found that were helpful:
1. https://r.qcbs.ca/workshop07/book-en/mixed-model-protocol.html
2. https://r.qcbs.ca/workshop07/book-en/step-3.-model-validation.html
3. https://r.qcbs.ca/workshop07/book-en/step-4.-interpretation-and-visualization.html
4. https://ourcodingclub.github.io/tutorials/mixed-models/
5. https://rpubs.com/daharo_calpoly/502695

Some important takeaways:
- "If your random effects are there to deal with pseudoreplication, then it doesn’t really matter whether they are “significant” or not: they are part of your design and have to be included."
- When doing model selection and comparing them, you cannot compare lmer models with lm models. Also, you need to make sure any of your mixed effects models are run with "REML=FALSE" so that the AIC can be correctly compared. But, once your final model is selected and it has fixed and random effects, you need to switch it back to "REML=TRUE".
- For post-hoc comparisons, 
