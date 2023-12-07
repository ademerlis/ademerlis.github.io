---
layout: post
title: Revisiting Pegasus (UM Supercomputer) after 4 years
date: '2023-05-24'
categories: coding
tags: [coding]
---

What has changed with Pegasus since 2019? Time to find out.

The basics: https://acs-docs.readthedocs.io/

also I found this article helpful: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8341507/

I did reset my password so I should be able to log in using that.

```{bash}
#log onto Pegasus
ssh and128@pegasus.ccs.miami.edu

#see what your bash profile looks like because its been four years
nano ~/.bash_profile

#PATH=$PATH:$HOME/bin
#export PATH
#export PATH=$PATH:~/software/local/parallel/bin
#alias and="cd /projects/scratch/transcriptomics/allysondemerlis"
#alias compute="bsub -P transcriptomics -Is bash"
```
Also I have used 140 GB of my 250 GB quota in my "nethome" user project space. 

```{bash}
#list files in home directory with sizes
ls -lh 

#total 2.9G
#drwxr-xr-x 2 and128 ccsuser          512 Sep 30  2019 logs
#drwxr-xr-x 4 and128 ccsuser         4.0K Mar  4  2021 scripts
#-rw-r--r-- 1 and128 ccsuser         2.4G Mar  4  2021 scripts.tar.gz
#drwxr-sr-x 7 and128 transcriptomics  512 Nov 29  2019 sequences
#drwxr-xr-x 4 and128 ccsuser          512 Sep 23  2019 software
#-rw-r--r-- 1 and128 ccsuser         445M Mar  4  2021 wound_healing.tar.gz
```

I don't remember making these .tar.gz files. I can't find any info on whether gunzipping them is resource-intensive and therefores should not be done on the login node. I guess to be safe I'll switch nodes (How do i do that again...)

After some googling, I don't think i necessarily have to "switch" to a computing node, it's just that when I submit jobs using bsub I specify the computing node or project scratch space or whatever it is. 

BUT - I just tried to move into my transcriptomics project space and it said "access denied." LOL. so first i need to get access again. But I think now they charge the PI of the project space for each user and amount of time used in the project space. I got an Early Career Researcher Award that gives me my own project space, so I just applied for that and hopefully I'll hear back soon. 
