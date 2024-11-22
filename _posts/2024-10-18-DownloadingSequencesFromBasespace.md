---
layout: post
title: Downloading sequences from Basespace
date: '2024-10-18'
categories: Analysis, Processing
tags: [coding]
---

NOTE: There is a great tutorial and overview from Danielle Becker in Dr. Hollie Putnam's lab on their GitHub lab notebook [here](https://github.com/Putnam-Lab/Lab_Management/blob/master/Bioinformatics_%26_Coding/Data_Mangament/SRA-Upload_Protocol.md).

I got these errors for some of the sequence files when I tried to upload them to NCBI:

<img width="844" alt="Screen Shot 2024-10-18 at 10 20 17 AM" src="https://github.com/user-attachments/assets/7c7d0ba2-cf22-48f5-9213-957fcc5ac44c">

And when I emailed NCBI, they said this:

<img width="624" alt="Screen Shot 2024-10-18 at 10 20 52 AM" src="https://github.com/user-attachments/assets/165df4c8-9250-4c65-9307-3e7d586275f0">

Michael Studivan also ran into this issue when he was uploading files to NCBI for a project. He said what solved it was finding the original files and uploading them directly.  

I tried that already from downloading them from the Basespace website to my local drive, then uploading them again. That didn't help. But what did seem to do something was if i gunzipped the file, then changed the metadata table so the file name was ".fastq" instead of ".fastq.gz", and then uploaded the file. However, now it looks like there is another file that doesn't work anymore.

<img width="815" alt="Screen Shot 2024-10-18 at 12 20 48 PM" src="https://github.com/user-attachments/assets/b20affa1-a8ff-4bec-a49c-794b8c14d244">

What's annoying is it doesn't seem to show me the full list of samples. So I really can't tell if it's all the Pcli samples, or just the ones listed there. 

I'm going to try to download all the Pcli files directly from Basespace using Michael's command line code ([here](https://github.com/mstudiva/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt)).

I think Michael's code is outdated, so I found the BaseSpace CLI tutorial on Illumina [here](https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview).

```{bash}
# if you have not previously, download BaseSpaceCLI
wget "https://api.bintray.com/content/basespace/BaseSpaceCLI-EarlyAccess-BIN/latest/\$latest/amd64-linux/bs?bt_package=latest" -O $HOME/bin/bs # this didn't work
chmod +x ~/bin/bs

#download via homebrew
brew tap basespace/basespace && brew install bs-cli

export PATH="/opt/homebrew/bin:$PATH"

echo '#!/bin/bash' > downloadReads.sh
echo 'bs download project --concurrency=high -q -n JA23031 -o . | grep "Pcli" ' >> downloadReads.sh
# -n is the BaseSpace project name and -o is the output directory

chmod +x downloadReads.sh

bs --version
# BaseSpaceCLI 1.6.1 -- built on 2024-07-23 at 09:33

bs auth
# has you go to the illumina website and sign in

./downloadReads.sh

```

the grep Pcli part didn't work but all the files are downloading very slowly. It looks like you can't specify which files to download 

```{bash}
bs list datasets --filter-term=".*Pcli.*"
#this lists only the basespace files for Pcli

bs download project -n JA23031 --exclude '.*Acer.*'
```

I also don't have space on my computer to download them locally. I have about 10 GB left of space on my computer. I can't decide whether it is worth it to try gunzipping and reuploading all the Pcli ones as fastq files to NCBI, or whether it is better to just try to redownload all the sample files from Illumina Basespace using the command line interface (CLI). I'm going to try to redownload all the files directly to an external hard drive and see how that goes.

Update: I think I successfully downloaded all 96 sequence files (.fastq.gz) onto the 2TB external hard drive on Nov 5.

