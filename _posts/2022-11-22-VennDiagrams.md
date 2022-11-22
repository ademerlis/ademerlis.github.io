---
layout: post
title: Venn Diagrams
date: '2022-11-22'
categories: Coding
tags: [Coding]
---

Today's goal is to make a venn diagram showing overlap of any significant DGEs between the time points. Similar to something like this from [Traylor-Knowles et al. 2021](https://www.frontiersin.org/articles/10.3389/fmars.2021.681563/full) and [Connelly et al. 2020](https://www.sciencedirect.com/science/article/pii/S0145305X20300094)

<img width="495" alt="Screen Shot 2022-11-22 at 9 12 45 AM" src="https://user-images.githubusercontent.com/56000927/203335783-8ce836fc-70a2-468f-99de-8e112d4881c9.png">
<img width="243" alt="Screen Shot 2022-11-22 at 9 13 45 AM" src="https://user-images.githubusercontent.com/56000927/203336001-245fa7d1-55f2-4a9f-92ce-0d770c6046ac.png">

Code from [Nick Kron](https://github.com/ademerlis/sctld_transcriptomics_2021/blob/main/orthofinder_mapping_NKron.Rmd) and [Mike Connelly](https://github.com/michaeltconnelly/EAPSI_Pocillopora_LPS/blob/master/Rmd/LPS_DESeq-venn_Pdam.Rmd)

I got this to work pretty easily (you have to manually calculate all the intersects which is kinda annoying... what's the point of the Venn Diagram package then)

```{r}
# need to manually calculate the totals and the intersects of each interaction (0-1, 0-2, 0-4, 1-2, 1-4, 2-4)

# use data frames anno_h0, anno_h1, anno_h2, anno_h4
# all DEGs is annotated_DEG_list

length(annotated_DEG_list$Row.names) #142 total

##Pairwise
h0_h1 <- intersect(anno_h0$Row.names, anno_h1$Row.names)
h0_h2 <- intersect(anno_h0$Row.names, anno_h2$Row.names)
h0_h4 <- intersect(anno_h0$Row.names, anno_h4$Row.names)
h1_h2 <- intersect(anno_h1$Row.names, anno_h2$Row.names)
h1_h4 <- intersect(anno_h1$Row.names, anno_h4$Row.names)
h2_h4 <- intersect(anno_h2$Row.names, anno_h4$Row.names)

##Triple
h0_h1_h2 <- intersect(h0_h1, anno_h2$Row.names)
h0_h2_h4 <- intersect(h0_h2, anno_h4$Row.names)
h1_h2_h4 <- intersect(h1_h2, anno_h4$Row.names)
h0_h1_h4 <- intersect(h0_h1, anno_h4$Row.names)
##Quadruple
h0_h1_h2_h4 <- intersect(h0_h1, h2_h4)

DGE_hour_venn <- draw.quad.venn(area1 = length(anno_h0$Row.names),
               area2 = length(anno_h1$Row.names),
               area3 = length(anno_h2$Row.names),
               area4 = length(anno_h4$Row.names),
               n12 = length(h0_h1),
               n13 = length(h0_h2),
               n14 = length(h0_h4),
               n23 = length(h1_h2),
               n24 = length(h1_h4),
               n34 = length(h2_h4),
               n123 = length(h0_h1_h2),
               n124 = length(h0_h1_h4),
               n134 = length(h0_h2_h4),
               n234 = length(h1_h2_h4),
               n1234 = length(h0_h1_h2_h4),
               category = c("Hour 0", "Hour 1", "Hour 2", "Hour 4"),
               fill = c("olivedrab3", "skyblue", "springgreen", "deepskyblue"),
               lwd = rep(1,4),
               cat.cex = rep(1.4, 4),
               #cat.fontfamily = rep("Arial", 4),
               cat.fontface = rep("plain", 4),
               alpha = c(0.3,0.3,0.3,0.3),
               cex = c(rep(1.4,5), 3, rep(1.4, 9)),
               fontface = c(rep("plain", 5), "bold", rep("plain", 9)),
               #fontfamily = c(rep("Arial", 15)),
               label.col = c(rep("black", 5), "red", rep("black", 9)),
               print.mode = "raw",
               sigdigs = 2,
               scaled = TRUE
) 
pdf("~/OneDrive - University of Miami/Cnidimmunity Lab/Wound_Healing_Project/2022 updates with Sami Beasley/manuscript figures and tables/VennDiagram_DGEallhours")
grid.draw(DGE_hour_venn)
```
<img width="499" alt="Screen Shot 2022-11-22 at 9 42 58 AM" src="https://user-images.githubusercontent.com/56000927/203342857-9c585a78-c04b-4fdd-b551-8a2a52168561.png">

When I look up that one gene that's shared among all time points, it is pdam_00003271 which was unannotated in the .gff file.

