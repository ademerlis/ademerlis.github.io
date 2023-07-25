---
layout: post
title: using Pcli transcriptome
date: '2023-07-25'
categories: coding
tags: [coding, Ch2_tempvariability]
---

So while the stringtie assembly is running for Acer, I am going to try to figure out how to use the Pcli file I downloaded (see [previous blog post](https://github.com/ademerlis/ademerlis.github.io/blob/master/_posts/2023-07-21_Pclitranscriptome.md)). 

First i gunzipped it so it would be .fna format. 

Then I just did "head" to see what the file looked like:

<img width="690" alt="Screen Shot 2023-07-25 at 9 59 24 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/36eafd47-8ef3-42b8-942b-337f2220eb55">

the "len=" and "path=[...]" are not normally part of a fasta file. 

So I'm confused about many things related to this -- so I have this transcriptome that was assembled de novo from a metatranscriptome by Avila-Magana et al. 2021. But we just have the .fna file, no gff file or anything. In reading their methods, it sounds like they used Kallisto to do pseudoalignment to a de novo reference transcriptome (so using their own sequences and then tring to align to that?) and then they used GO_MWU to get gene ontology terms. But how did they go from the transcriptome to GO_MWU, because that requires gene IDs with GO terms. 

Ok in their methods they say this: "Amino acid sequences were predicted from the Coral host and Symbiodiniaceae transcriptomes by using Transdecoder (v2.0.1)87. Orthologous groups of protein sequences amongst the three coral and their associated three photosymbiont species were determined with the OrthoFinder (v2.4.0)16 bioinformatics tool, using default parameters. Using reciprocal best-hits via BLAST all-v-all algorithm, Orthofinder determined the number of conserved putative orthologues among the three coral and three photosymbiont species. By a custom python script, the ortholog differential gene expression (logFC) was retrieved for each expression dataset per species against the other two species."

That sounds like how they did their "annotations". Because what I'm missing is basically the connection from the host transcriptome to any sort of identifiable information. 

Or what is we do something like this:

" To generate gene ontology annotations from our 169 assemblies we utilized the annotation software EggNog mapper (Huerta-Cepas et al., 2017, 170 2019) through the online web portal against the eukaryote database. To annotate individual 171 contigs we used blastp against the swiss-prot uniprot database (The UniProt Consortium et al., 172 2021)." 
(from [Dimos et al. 2022](https://doi.org/10.1101/2021.04.28.441826)).

Is that how you make a gff file? 
