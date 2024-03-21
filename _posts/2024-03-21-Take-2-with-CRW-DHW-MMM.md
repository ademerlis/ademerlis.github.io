---
layout: post
title: Take 2 with CRW DHW MMM
date: '2024-03-21'
categories: Analysis, Processing
tags: [Ch4_AcerCCC, coding]
---

I am revisiting this post from last year, where I was trying to download data from NOAA's CRW to calculate MMM and DHW for my reef sites. [Post 
here](https://github.com/ademerlis/ademerlis.github.io/blob/main/_posts/2023-12-04-DownloadingSSTCoralWatchData.md)

I found Dr. Ana Palacio's code on [her GitHub](https://github.com/anampc/DHW_Uva/tree/master) that documents really well her scripts for downloading the data from NOAA's database. 

First, I downloaded the .nc file here: https://www.star.nesdis.noaa.gov/pub/sod/mecb/crw/data/5km/v3.1_op/climatology/nc/

Then, I needed to download a package for the terminal that allows me to read .nc files.

```{bash}
brew update
brew install netcdf

ncdump -h ct5km_climatology_v3.1.nc
```

I can see in the header that I think it includes all lat, lon coordinates. 

Now, I read it into R using Ana's code:

```{r}

```
