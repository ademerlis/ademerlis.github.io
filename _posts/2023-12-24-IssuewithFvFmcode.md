---
layout: post
title: Issue with FvFm code
date: '2023-12-24'
categories: Coding, IPAM
tags: [Coding, IPAM]
---

I was looking at the full dataset of all the IPAM data to see if there were any issues, because when I started trying to remove NAs from the GLM for FvFm treatment data, I started getting errors. 

I noticed when I was looking at the different time points that six Acer are missing initial FvFm data, and when I look at the csv files that I exported from the IPAM, they have values that I guess weren't imported into R using the [IPAM2R code](https://github.com/ademerlis/ademerlis.github.io/blob/main/_posts/2023-01-21-fixingIPAM2Rcode.md). I did edit the code because it wasn't working, maybe that's why it is importing things in correctly? 

Ok, upon further inspection, the Acer that are missing the initial time point fvfm are all from one photo: photo number 4 ("4.csv") from 2022-03-16. All the Acer that are missing the final time point fvfm are from tank 7 (2022-04-20). For Pclivosa there are several final time point ones that are missing (2022-04-20). Some of them are actually blank in the original IPAM csv files, but some of them have values that should've been read in.

So what happened there? I guess let me screenshot to keep record of which ones were missing, and then re-import everything using the IPAM2R code and see if it corrects itself? 

![Screen Shot 2023-12-24 at 8 17 01 AM](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/d4da8990-c05a-41dd-9b2c-c826ce4cf5df)

UPDATE: I figured out the issue. The IPAM2R code has a line in it that uses the reshape2 function, and some of my csv files had two rows of identical IPAM data instead of one row - one row is what the IPAM2R code expects. I guess if you export your IPAM data twice on accident it will add another row to your csv file. But, when the function in R tries to reshape the data frame after tidying the AOI and other variable info, it basically extracts the wrong information and the Y(II) values end up being incorrect. For the Acer samples that were messed up, they had valuess of 2.00 instead of 0.646 etc that they were supposed to be. So I added a if-then statement within the IPAM2R data which basically removes the second row if there are two rows present in the csv file. 

When I looked through the updated dataset that was created in R, I no longer had so many NAs as I did before. I cross-validated the values from the original csv files to make sure they matched, and they did. Yay problem solved!
