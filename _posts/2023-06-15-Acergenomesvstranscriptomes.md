---
layout: post
title: Acer genomes vs transcriptomes
date: '2023-06-15'
categories: coding
tags: [coding, CCC_ch4, temperaturevariability2023]
---

I wanted to make a separate post dedicated to Acer because it seems particularly confusing. 

First, a lot of publications have been using the Libro et al. 2013 Acropora cervicornis transcriptome (https://www.dropbox.com/s/wovxoi2dxcz9kv2/a.cervicornis_2014.zip?dl=0&file_subpath=%2Fa.cervicornis), which it is important to note that this is a holobiont transcriptome, so it contains everything the Acer hosts as well. This paragraph from the supplemenatry information of Parkinson et al. 2018 describes it here: 


![image](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/19aa9b58-24cf-4d30-8d7a-f80d57ec0337)


In Lesneski et al. 2022 (in prep: https://www.biorxiv.org/content/10.1101/2022.03.29.486305v1), the authors state that the Libro et al. 2013 paper was based on RNA sequence data from both Acer and Apal, and this risks generating "mosaic contigs". And, the transcriptome has a relatively low "N50, suggesting a paucity of full-length transcripts." Lesneski et al. generated two "better" transcriptome assemblies for Acer from Belize. Fasta files located here: https://osf.io/u9nq8/

Baums and Kitchen have an Acer genome hosted on Galaxy (not technically published but publicly available?) 
https://usegalaxy.org/published/page?id=2f8d21c73f8501e2
https://usegalaxy.org/published/history?id=1f8678b27ae56467


Ok, so a common theme keeps coming up. In Libro et al, Lesneski et al, and Baums et al, there are no gtf annotation files. Am i supposed to build that? 

Looking at this STAR tutorial (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4631051/), this is what they said regarding transcript/gene annotations in GTF format:

"The Basic Protocol uses transcript/gene annotations in GTF format. The gene annotations allow STAR to identify and correctly map spliced alignments across known splice junctions. While it is possible to run the mapping jobs without annotations, it is not recommended. When gene annotations are not available, use the 2-pass mapping described in Alternate Protocol 2."

I think what is confusing me too is looking at Natalia's index script (https://github.com/China2302/SCTLD_RRC/blob/main/hpc/STAR_index.sh), she has a file for Ofav genome that is "genomic.fna" as the fasta file, and then "genomic.gtf" as the GTF file. But, when I look at the files available for the Acer genome from Iliana Baums, there are different .fa files (CDS, protein) versus the scaffold .fasta file. Which one do I use to build the index? Why isn't there a gtf file? there is a text annotation file, but can I convert that to gtf?

Alternatively, should I try the Matz / Bowtie2 route? That seems to be another popular pipeline that publications use for Tag-seq based transcriptomes/genomes in particular. 

Update: I found these two articles which say that gff3 and gtf are similar files and can be used in STAR under the --sjdbGTFfile flag. However, one source suggests that converting the gff to gtf using "gffread" is recommended (idk why but I think it has to do with making sure the information is presented correctly for STAR to use it). https://biohpc.cornell.edu/doc/RNA-Seq-2018-Lecture1.pdf and https://goenomics.com/glossary.html

See https://github.com/ademerlis/ademerlis.github.io/blob/master/_posts/2023-06-21_STARalignmenttroubleshooting.md for more details
