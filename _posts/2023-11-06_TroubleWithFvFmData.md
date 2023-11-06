---
layout: post
title: Trouble with Fv/Fm Data
date: '2023-11-06'
categories: Coding
tags: [Coding]
---

I'm revisiting the IPAM data from the Chapter 2 stress-hardening experiment, and I'm having a really hard time wrapping my head around it.

I think I've finally come up with a workflow for using GLMMs to test significance of percentage/proportion (0<x<1 bound data): [see here](https://github.com/ademerlis/temperaturevariability2023/blob/main/physiology/photosynthetic_efficiency/fvfm.Rmd). I basically follow Kevin Wong's GLM code. 

**My major confusions right now are:**
1. When I calculate normalized CBASS Fv/Fm values based on the final Fv/Fm measurement at the end of the treatment period, I get some negative numbers and numbers > 1. There are too many of them to filter out, so I don't know how to statistically test these.
2. I want to test the significance of the dose-response curves (DRC) for the CBASS values, but instead of ED50 I want to test the slopes of the curves to see if those are significant. In the drc model code for R, you can specify parameters (i.e. ll.3 = "hill", "max", "ED50") and I believe "hill" is the "hill slope". I can extract those values from the model, but idk what they mean and how to statistically test them.

