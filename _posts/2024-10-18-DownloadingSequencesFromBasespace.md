---
layout: post
title: Downloading sequences from Basespace
date: '2024-10-18'
categories: Analysis, Processing
tags: [coding]
---

I got these errors for some of the sequence files when I tried to upload them to NCBI:

<img width="844" alt="Screen Shot 2024-10-18 at 10 20 17 AM" src="https://github.com/user-attachments/assets/7c7d0ba2-cf22-48f5-9213-957fcc5ac44c">

And when I emailed NCBI, they said this:

<img width="624" alt="Screen Shot 2024-10-18 at 10 20 52 AM" src="https://github.com/user-attachments/assets/165df4c8-9250-4c65-9307-3e7d586275f0">

Michael Studivan also ran into this issue when he was uploading files to NCBI for a project. He said what solved it was finding the original files and uploading them directly.  

I tried that already from downloading them from the Basespace website to my local drive, then uploading them again. That didn't help. But what did seem to do something was if i gunzipped the file, then changed the metadata table so the file name was ".fastq" instead of ".fastq.gz", and then uploaded the file. However, now it looks like there is another file that doesn't work anymore.

<img width="815" alt="Screen Shot 2024-10-18 at 12 20 48 PM" src="https://github.com/user-attachments/assets/b20affa1-a8ff-4bec-a49c-794b8c14d244">

What's annoying is it doesn't seem to show me the full list of samples. So I really can't tell if it's all the Pcli samples, or just the ones listed there. 
