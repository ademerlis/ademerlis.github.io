---
layout: post
title: TopGO = what does it all mean?????
date: '2023-02-08'
categories: RNA Sequencing
tags: [RNA, Sequencing]
---

I am not sure how to interpret the results of the topGO analysis.

For example, in the results table for upregulated gene GO terms for hour 1:

![Screen Shot 2023-03-22 at 2 33 20 PM](https://user-images.githubusercontent.com/56000927/227004112-1b354047-b24b-46a2-9d3c-3c6db7625ca4.png)

What does the annotated column mean? There are only 5 significant DGEs that are being inputted into this GO enrichment. When you look at the significant column, the maximum number is 1, which means that only 1 significant DGE is part of that GO term (I think). So is the annotated column the total number of genes in the gene universe that are annotated to that GO term? Is that useful or relevant? Maybe that's what's used to calculate the expected probability of that GO term being represented then? 

Also, this figure I made: 
<img width="709" alt="Screen Shot 2023-03-22 at 2 37 55 PM" src="https://user-images.githubusercontent.com/56000927/227004705-8afaabb9-df08-4a69-b85a-70cdd6d73a2c.png">

Are annotated genes even relevant in this figure? that has nothing to do with my dataset (I don't think anyways - other than the gene universe was built from the total reads from my samples). In Traylor-Knowles et al. 2021, we used "Number of Protein Identifiers." Does that mean number of annotated genes within each GO term? 

I am so confused what this all means.

I found this code from avrilomics which helps find the list of genes "annotated" (that's the word avrilomics used, but idk if that's the appropriate use of the word given the context) to a particular set of GO terms (http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html):

```{r}
myterms = c("GO:0007610", "GO:0014070", "GO:0045910")
mygenes <- genesInTerm(myGOdata, myterms)
for (i in 1:length(myterms))
   {
       myterm <- myterms[i]
       mygenesforterm <- mygenes[myterm][[1]]
       mygenesforterm <- paste(mygenesforterm, collapse=',')
       print(paste("Term",myterm,"genes:",mygenesforterm))
     }
```

Ok, from this source I think I figured it out (https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html):

- Annotated: number of genes (in our gene list) that are annotated with the term
- Significant: Number of significantly DE genes annotated with that term (i.e. genes where geneList = 1)
- Expected: Under random chance, number of genes that would be expected to be significantly DE and annotated with that term
- raw.p.value: P-value from Fisher’s Exact Test, testing for association between significance and pathway membership

Also, based on the above source, the scoring I used for the GO analysis ("1" if padj < 0.05 and "0" if padj > 0.05) means that I need to use the Fisher Exact test, not the K-S test, because the K-S test compares the p-values themselves. 
