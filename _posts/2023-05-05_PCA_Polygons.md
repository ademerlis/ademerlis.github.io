---
layout: post
title: Creating Polygons for PCAs
date: '2023-02-08'
categories: RNA
tags: [RNA]
---

![Screen Shot 2023-05-05 at 11 37 11 AM](https://user-images.githubusercontent.com/56000927/236503510-e0fbf806-a960-4bc4-82f9-513c279b28e3.png)

Figure 1. Principal component analysis of wounding (orange triangle) versus control (blue circle) treatments at hour zero (A), one (B), two (C), and four (D) post-injury. Counts normalized using variance stabilizing transformation (VST) from DESeq2 (Love et al. 2014). 

I want to create polygons with these points so that it's easier to see the clustering. 

Ok that was surprisingly easy, just use geom_polygon() 

![Screen Shot 2023-05-05 at 11 57 14 AM](https://user-images.githubusercontent.com/56000927/236507879-544c798e-8e02-40de-8024-1735236dd8bf.png)

Maybe I can try to add the PERMANOVA results on there too? 

All the code for my PERMANOVAs says that none of the PCAs are significant, but hours 2 and 4 look like they would be significant.

After googling this for a couple hours, I'm going to just leave it as is and not add stats to the PCA because I don't think it's necessary. Also I am not sure if the tests I even did were correct.

Update May 12, 2023:

I found this code in my old notes from the SCTLD Jamboree, could be useful for calculating significant distances between polygons:

```{r}
# Calculate distances among samples
sampleDists <- dist(t(assay(dds_vst_ofav)), method = "manhattan")
sampleDistMatrix <- as.matrix(sampleDists)

# Calculate MDS
mds <- as.data.frame(colData(dds_vst_ofav)) %>% 
  cbind(cmdscale(sampleDistMatrix))

# Calculate MDS and use eigenvectors to determine proport
mds_eig <- cmdscale(sampleDistMatrix, eig = TRUE)
mds_eigenvectors <- data.frame(mds_eig$eig) %>% 
  mutate(prop_var = mds_eig.eig / sum(mds_eig.eig))

# Calculate Treatment centroids for plotting
mds_trmt <- mds %>%
  group_by(condition) %>%
  dplyr::summarise(c1 = mean(`1`), c2 = mean(`2`)) %>%    
  full_join(mds)

# Plot with spiders
ggplot(mds_trmt, aes(fill =condition)) +
  #stat_ellipse(aes(x = `1`, y = `2`, color = condition, fill = condition), geom = "polygon", type = "norm", alpha = 0.0) + 
  # sample-centroid spiders paths
  geom_segment(mapping = aes(x = `1`, y = `2`, xend = c1, yend = c2),
               lwd = 0.25, col = "dark grey") +
  # treatment centroid points
  geom_point(size = 3, aes(x = c1, y = c2, color = condition), fill = "black", shape = 21, stroke = 2, show.legend = TRUE) +
  # sample points
  geom_point(size = 3, aes(x = `1`, y = `2`, color = condition), stroke = 0.5, show.legend = FALSE) +
  scale_color_manual(values = (c("#7CBAF5", "#AD161A"))) +
  theme_bw() +
  labs(x = "MDS1", y = "MDS2")

```
