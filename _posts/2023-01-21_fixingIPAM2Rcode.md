---
layout: post
title: Figuring out why R code for importing PAM data stopped working
date: '2023-01-21'
categories: Coding
tags: [Coding]
---

I realized that there were some issues with the RStudio on my AOML desktop versus my laptop. First, I had updated the RStudio on my computer, so it messed up one of the packages that was needed for Ross Cunning's IPAMtoR custom program. I tried to install devtools to install this package (called "joeyr") but AOML wasn't letting me do that for some reason, or something else needed to be updated or something. I got it to work on my laptop though. Then, once I actually got it to work, I have been having issues importing the csv files that were downloaded from the IPAM software. 

The error I keep getting (just for two of the files which is weird) is "Error in if (is.na(dplyr::last(df))) df <- df[, -(ncol(df))] : 
  the condition has length > 1"
  
I thought maybe there was a download issue from OneDrive, because everything is stored on OneDrive. But even after downloading everything, still it is not working. So maybe it is something wrong with the import_ipam function. I am going to try to comment out the line "if (is.na(dplyr::last(df))) df <- df[,-(ncol(df))]", which is meant to get rid of trailing NAs, and see if that will let me import them.

Now the error I get is "Error in dim(ordered) <- ns : 
  dims [product 1] do not match the length of object [0]"

I can't figure out why this isn't working. Maybe I need to save the files as .csv files again?

<u> Jan 21st 2023 update </u>  

I installed an older version of RStudio (desktop from 7.1.2022) and tried to load the ipam files again, and march 16 and april 20 still didn't work (i get the same error from before). So I am going to download the files from google drive again and see if that helps in any way. maybe the files got reformatted somehow.

I did this for march 16 and it still didnt work...

I also tried to isolate the files that the error code was popping up for (it seemed like sequentially some of the files worked, but then there would be one that the code would get hung up on). But when testing csv files individually, they all worked and there were no issues.

Maybe I need to ask Ross, Carly, Rich, or Alex if they have ever come across this error.

OK WAIT I GOT IT TO WORK!!!!

something about the if-then statement, R was getting tripped up on. So i replaced the lines:
# Get rid of trailing NA
    if (is.na(dplyr::last(df))) df <- df[,-(ncol(df))]
    
with:
# Get rid of trailing NA
    df <- df %>% select_if(~ !any(is.na(.)))
    
And I can now import march 16 and april 20!!!!! WOOOO!
