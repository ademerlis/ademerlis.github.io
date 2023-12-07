---
layout: post
title: installing MultiQC on pegasus
date: '2023-06-01'
categories: coding
tags: [coding]
---

I want to install multiqc so I can look at all the fastQC results at once and compare them to one another. 

https://multiqc.info/#:~:text=Install%20from%20the%20Python%20Package,installation%20instructions%20for%20more%20help.

```{bash}
module load py-pip/20.2
pip install multiqc
```

this is taking awhile but it might work?

Ok I think it worked, but I got this warning for a bunch of the related packages that got installed along with multiqc:

WARNING: The script pygmentize is installed in '/nethome/and128/.local/bin' which is not on PATH.
  Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.

I think PATH is something related to my bash_profile? 

I'm just going to try running multiqc as is and see if it works. 

Ok navigating into the folder with all the fastqc results and running "multiqc ." did not work. It said "command not found."

I need to figure out how to load a local package.

I found this tutorial on a cutadapt tutorial website: https://gensoft.pasteur.fr/docs/cutadapt/1.18/installation.html

 " You can then run the program like this:
  ```{bash}
  ~/.local/bin/cutadapt --help
  ```
  If you want to avoid typing the full path, add the directory $HOME/.local/bin to your $PATH environment variable."
  
So I need to add the .local/bin directory to $PATH so that I can run it without having to use the whole path.

Ok, I did this and it worked to run multiqc:
```{bash}
nano ~/.bash_profile
export PATH=$PATH:/nethome/and128/.local/bin

#then save it and exit nano

source .bash_profile
#this runs the bash profile to update the paths
```

Then I navigated the the directory with all the fastqc results and ran "multiqc ."

I will transfer the multiqc_report.html to my local drive so I can open it and view it.

this code worked from transferring from pegasus to local (first navigated on local to the folder i wanted the report to go in):
```{bash}
Allysons-MacBook-Pro-2:AcerCCC allysondemerlis$ scp and128@pegasus.ccs.miami.edu:/scratch/projects/and_transcriptomics/Allyson_CCC/multiqc_report.html .
```




  
