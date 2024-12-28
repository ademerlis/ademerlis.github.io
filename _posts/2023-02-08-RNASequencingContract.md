---
layout: post
title: RNA Sequencing Contract with Dr. Michael Studivan
date: '2023-02-08'
categories: [RNA-Seq, Budgeting, Dissertation]
tags: [RNA, Sequencing, RNA-Seq, Stress-Hardening, Reciprocal Transplant, Michael, Budgeting]
---

Prepping metadata from RNA extractions for 2022 Stress-Hardening Experiment and 2022-2023 Urban Coral Reciprocal Transplant (Carly+Rich Experiment)

- NOAA Sequencing contract is with UT Austin GSAF:
https://wikis.utexas.edu/display/GSAF/Home+Page

- Need to prepare plate maps with Sample IDs and their locations on a 96-well plate. 
https://docs.google.com/spreadsheets/d/1Ed0iTtE-sqT7BKSp2zzZTqw-9wOGhohaJa2O7GcZdzU/edit#gid=0

- We also calculated dilutions from concentrations of RNA determined using Qubit.
Stress-Hardening project: https://docs.google.com/spreadsheets/d/1P1Lgnrm11YLqmGaR0SajlSMjOQAv6jyE/edit#gid=1845638551
Carly/Rich Reciprocal transplant: https://docs.google.com/spreadsheets/d/1N3xmNduMS7OGXOdKx_oDRB-9T7m19WpMLewaFyqhyfY/edit#gid=0

UT Austin prefers to receive 25 uL of sample, and Michael determined based on the quantities of RNA that we extracted that we could get away with aiming for 500 ng of RNA per sample. We eluted RNA in either ~20 or ~40 uL of nuclease-free water (depending on whether or not then went through the Zymo Clean and Concentrate-5 Kit). That's important to keep track of in a spreadsheet too for all the metadata.

As for the sequencing contract itself, we have 160 samples that we are submitting for sequencing. You have to specify sequencing type and read depth (coverage) as well. 

We chose library prep type = TagSeq
- For TagSeq coverage, we went for high coverage (7-10 M). 
- We want single-end 100 base pairs.

Cost per Sample from UT Austin website: 
- RNA-Seq Library Preparation, with no enrichment = $178.92
- RNA-Seq Library Preparation, with poly(A) Enrichment = $260.39
- RNA-Seg Library Preparation, with kibosomal kemov. = $380.62
- Small RNA Library Preparation = $178.92
- TagSeq High Coverage (7-10M) Library Preparation = $45.18
- TapSen Standard Coverage (3-5M) thrar, Preparation = $45.18

For the type of sequencer, we chose NovaSeq. MiSeq is not appropriate for this (long read but low coverage or low depth or something). NexSeq has a lot of problems with TagSeq.

For NovaSeq, there are different options within that as well. There are different flow cell types, which differ in overall number of reads that they generate. 
NovaSeq S1 PE50, SR100 = $8300 (internal UT Austin cost) 1.4 billion reads total
NovaSeq S2 PE50, SR100 = $11,141 (internal UT Austin cost) 3.7 billion reads total

We ended up doing the combination of TagSeq High Coverage on NovaSeq S1 --> totalled around $17,000 for 160 samples.
