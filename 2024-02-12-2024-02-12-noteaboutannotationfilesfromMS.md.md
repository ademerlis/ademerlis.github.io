---
layout: post
title: 2024-02-12-noteaboutannotationfilesfromMS.md
date: '2024-02-12'
categories: Analysis
tags: [Ch2, coding]
---

I realized that when Michael re-annotated the Locatelli unpublished Acer genome in December, he updated 
the names of everything (which no longer match my counts matrix). Basically when we first started this 
pipeline, all the Acropora gene names were "Acropora_" before the gene number. But in the recent 
re-annotation, Michael named everything "Acervicornis" with no underscore. I ran into issues when I was 
trying to run the GO_MWU scripts because I was using the "Acervicornis_iso2go.tab" file but trying to 
match the "Acropora_" gene names (which none of them matched of course). 

So I transformed the new iso2go.tab file so that everything was consistent and named like "Acropora_". 

Today I just realized that in the KOG MWU script and in the annotating DGEs scripts for Acer, I was still 
using the old files that Michael had made (I realized this because the iso2geneName.tab and 
iso2kogclass.tab names are Acropora_, not Acervicornis). 

Okay but now looking back at his old commit versions on the Acer annotated transcriptome repository on his 
GitHub, it looks like the naming convention for gene names used to be "Acropora1921" for example, which 
means it didn't have an underscore and it also didn't have zeros before the numbers (which he then added 
when he changed it to Acervicornis).

Just to be safe, I'm going to redownload all the annotation files and rerun everything that uses those 
files so that way I can be completely sure its all updated.
