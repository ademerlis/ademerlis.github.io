---
layout: post
title: Uploading files to Pegasus
date: '2023-05-31'
categories: [Bioinformatics]
tags: [Pegasus, bash, CCC]
---

After trying several different codes, this is what I got to work:

```{bash}
Allysons-MacBook-Pro-2:2TB allysondemerlis$ scp -r Allyson_CCC/ and128@pegasus.ccs.miami.edu:/scratch/projects/and_transcriptomics/
```

First, you need to navigate to the location of the files you want to transfer to pegasus. In this case it's the "2TB" external hard drive. You have to go up two file directories out of Users to get there.

Then you specify the project space within your pegasus account, and then it will ask you to enter the password. 

And that's it!

I first did the Coral City Camera samples (N=20) and that took about 30 min. 

Now I'll move over the stress-hardening ones and then work on the code to do FastQC.
