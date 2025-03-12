---
layout: post
title: Downloading sequences from Box
date: '2025-03-07'
categories: Bioinformatics
tags: [Pegasus]
---

I tried some different methods for downloading the Chapter 3 Reciprocal Transplant experiment .fastq.gz sequences from Box, and realized that there isn't a way to download/upload the raw sequences from Box without locally downloading them to my computer (which would take 100s of GBs of space). So, I think my best bet is to use the computer in the NTK lab to download them locally, then upload them to Pegasus.

Things that I tried that didn't work: downloading some sort of Box CLI and trying to remotely access it via Pegasus. 

```{r}
curl -L https://downloads.boxcdn.net/boxcli/latest/BoxcliInstaller.run | bash\n
brew install rclone\n
nano ~/.bash_profile
source ~/.bash_profile
ls
brew install rclone\n
rclone config
rclone lsd
rclone lsd Box:/
less ~/.config/rclone/rclone.conf
scp ~/.config/rclone/rclone.conf and128@pegasus.ccs.miami.edu:/scratch/projects/and_transcriptomics/ #tried moving the configuration file from my computer to pegasus

# on Pegasus
module load rclone #rclone is already installed, but it doesn't have Box available
module list
rclone -V
rclone config
mv rclone.conf ~/.config/rclone/
rclone lsd box:
touch ~/.config/rclone/rclone.conf
```

There is a [section about remote transferring on the Pegasus website](https://acs-docs.readthedocs.io/services/2-transfer.html). 

It also still says you can use FileZilla, but when I try going to the FileZilla website on the "UMiamiWireless" wifi on campus, the page is blocked. I need to be on this wifi to access Pegasus, so I don't know a workaround for that.

(wait can't I just download FileZilla locally when on a different wifi and then try it)...

Regardless, I think using FileZilla would still require locally downloaded files from Box.

## March 12 Update:

I successfully downloaded sequences from Box drive to 2TB hard drive using the Cnidarian Immunity mini mac. It took about half a day (once I plugged in the ethernet cable).

I used "rsync" to have it give me updates in the terminal for each file, which i found really helpful.

this worked:
```{bash}
rsync -ah --progress --partial . "/Volumes/2TB/Ch3_RT"
```

Today, I'm using rsync again to transfer all files from 2TB to Pegasus scratch space (and_transcriptomics).

this is working:
```{bash}
rsync -ah --progress --partial -e ssh "/Volumes/2TB/Ch3_RT" and128@pegasus.ccs.miami.edu:/scratch/projects/and_transcriptomics/reciprocaltransplant/raw_seq_files
```
![Screenshot 2025-03-12 at 2 34 36 PM](https://github.com/user-attachments/assets/be313172-dda8-4a7f-9690-a7eb8b15569b)

What's weird is the .gz files on Pegasus are now green instead of red. But when i check the md5sum values for the original file and the transferred file, they match. So I guess there's no loss in file integrity, idk what else it could be.

![Screenshot 2025-03-12 at 2 35 06 PM](https://github.com/user-attachments/assets/5a43ac23-53b2-4466-aa1d-b4070c1fd0d2)


When I run this, I get this: "gzip compressed data, extra field"
```{bash}
file Pstr-Dec2022-164_S46_L001_R1_001.fastq.gz
```
So whatever, it's fine I guess.
