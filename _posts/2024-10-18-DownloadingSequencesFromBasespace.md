---
layout: post
title: Downloading sequences from Basespace
date: '2024-10-18'
categories: Bioinformatics
tags: [Basespace, NCBI, SRA]
---

NOTE: There is a great tutorial and overview from Danielle Becker in Dr. Hollie Putnam's lab on their GitHub lab notebook [here](https://github.com/Putnam-Lab/Lab_Management/blob/master/Bioinformatics_%26_Coding/Data_Mangament/SRA-Upload_Protocol.md).

I got these errors for some of the sequence files when I tried to upload them to NCBI:

<img width="500" alt="Screen Shot 2024-10-18 at 10 20 17 AM" src="https://github.com/user-attachments/assets/7c7d0ba2-cf22-48f5-9213-957fcc5ac44c">

And when I emailed NCBI, they said this:

<img width="500" alt="Screen Shot 2024-10-18 at 10 20 52 AM" src="https://github.com/user-attachments/assets/165df4c8-9250-4c65-9307-3e7d586275f0">

Michael Studivan also ran into this issue when he was uploading files to NCBI for a project. He said what solved it was finding the original files and uploading them directly.  

I tried that already from downloading them from the Basespace website to my local drive, then uploading them again. That didn't help. But what did seem to do something was if i gunzipped the file, then changed the metadata table so the file name was `.fastq` instead of `.fastq.gz`, and then uploaded the file. However, now it looks like there is another file that doesn't work anymore.

<img width="500" alt="Screen Shot 2024-10-18 at 12 20 48 PM" src="https://github.com/user-attachments/assets/b20affa1-a8ff-4bec-a49c-794b8c14d244">

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

The grep Pcli part didn't work but all the files are downloading very slowly. It looks like you can't specify which files to download 

```{bash}
bs list datasets --filter-term=".*Pcli.*"
#this lists only the basespace files for Pcli

bs download project -n JA23031 --exclude '.*Acer.*'
```

I also don't have space on my computer to download them locally. I have about 10 GB left of space on my computer. I can't decide whether it is worth it to try gunzipping and reuploading all the Pcli ones as fastq files to NCBI, or whether it is better to just try to redownload all the sample files from Illumina Basespace using the command line interface (CLI). I'm going to try to redownload all the files directly to an external hard drive and see how that goes.


**Update:** I think I successfully downloaded all 96 sequence files (.fastq.gz) onto the 2TB external hard drive on Nov 5.
When they download from Basespace, each sequence file is in its own subdirectory. For uploading to NCBI, all the sequence files have to be in one directory, no subfolders. So, I ran these two lines of code from Michael to easily fix that:

```{bash}
find . -name '*.gz' -exec mv {} . \;
rmdir SA*
```

## Nov 25:
I finally finished uploading all the sequence files using ftp from the external hard drive, and I am still getting the same error on NCBI.

<img width="500" alt="Screen Shot 2024-11-25 at 11 33 46 AM" src="https://github.com/user-attachments/assets/7a430a66-6dc9-4c5a-9e20-b9fcc9211838">

I have no idea what else to do besides gunzipping all the files and then re-uploading them as .fastq files instead. Maybe I should just try this for the Pcli first and see if it accepts the Acer as is.

## December 2

I started going through and running "gunzip" to each of the files that NCBI listed as an error in the screenshot above. It had me do this twice for a number of Pcli samples, and after re-uploading those as .fastq files, I got this list of flags (which I was anticipating at this point): (note that this is just a screenshot of the first few listed -- all of the samples are listed on the NCBI website as having this same error)

<img width="500" alt="Screen Shot 2024-12-02 at 2 17 09 PM" src="https://github.com/user-attachments/assets/f9f7455d-682b-4799-81e8-ab1c57498a5f">

aka all the files are "corrupted." 

Which I don't understand since I downloaded them straight from Basespace. So now I'm running "gunzip" on all of the files, and seeing if I re-upload them all as .fastq files, whether they will be flagged. It looks like based on the recent list that the .fastq files for the PCli samples I recently did were successfully uploaded without an error. But we shall see once I do all the rest of them.



