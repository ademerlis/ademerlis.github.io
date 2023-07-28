---
layout: post
title: Pcli Salmon transcriptome
date: '2023-07-26'
categories: coding
tags: [coding, Ch2_tempvariability]
---

So I installed Salmon locally onto my scratch space by downloading the Linux binary package:

```{bash}
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz
tar -xzvf salmon-1.5.2_linux_x86_64.tar.gz
cd salmon-1.5.2_linux_x86_64/bin/
#salmon is green so that means it is already compiled
#add to PATH
nano ~/.bash_profile
source ~/.bash_profile
#check it works:
salmon-1.5.2_linux_x86_64/bin/salmon --verson
```

This salmon code worked:
```{bash}
#!/bin/bash
#BSUB -J Pcli_transcriptome_salmon_index
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -o Pcli_transcriptome_salmon_index%J.out
#BSUB -e Pcli_transcriptome_salmon_index%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

${and}/programs/salmon-1.5.2_linux_x86_64/bin/salmon index -t ${and}/genomes/Pcli/clean_Pcli_transcriptome_final.fasta -i ${and}/genomes/Pcli/Pcli_transcriptome_index
```

Output:
![Screen Shot 2023-07-26 at 9 17 51 AM](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/1acc45ef-e328-483d-a193-6ba5535cc124)

Following [this Salmon tutorial](https://combine-lab.github.io/salmon/getting_started/):

```{bash}
#!/bin/bash
#BSUB -J Pcli_transcriptome_salmon_quant
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -o Pcli_transcriptome_salmon_quant%J.out
#BSUB -e Pcli_transcriptome_salmon_quant%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/AS_pipeline/3_trimmed_fastq_files/Pcli_fastq_files"

data=($(ls *.gz))

for sample in ${data[@]} ;

do \
${and}/programs/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${and}/genomes/Pcli/Pcli_transcriptome_index -l A -1 ${sample} -p 8 --validateMappings -o ${and}/Ch2_temperaturevariability2023/AS_pipeline/4_Pcli_specific/salmon_quant_files ; \
done
```

I got this error: "/scratch/projects/and_transcriptomics/programs/salmon-1.5.2_linux_x86_64/bin/salmon quant was invoked improperly"

I searched this on ChatGPT (Exception : [std::bad_alloc] was also part of the error message) and it said that it could be a memory allocation issue. So i changed the script to bigmem with 16 nodes, and that still got the same error message.

I tried to run STAR index instead and see if that worked but I kept getting a fatal error due to not enough RAM. Salmon is supposed to be more efficient than STAR anyways, but I don't understand why the Salmon quant keeps failing. 

STAR error:
<img width="654" alt="Screen Shot 2023-07-27 at 10 30 38 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/2c0e1eba-928a-45ff-abfe-450aa557abc5">

How much RAM allocation is allowed on Pegasus?

```{bash}
[and128@login4 scripts]$ free -h
              total        used        free      shared  buff/cache   available
Mem:           125G         18G         99G        254M        7.2G        105G
Swap:          4.0G        3.2G        779M
```
Mem is the amount of physical memory (RAM) available, and swap is the amount of virtual memory for when RAM is exhausted.

So total is... ~ 100 GB. 

For STAR, the error said: EXITING because of FATAL PARAMETER ERROR: limitGenomeGenerateRAM=31000000000 is too small for your genome
SOLUTION: please specify --limitGenomeGenerateRAM not less than 511308027488 and make that much RAM available.

It is requesting 511308027488 bytes, which = ~500 GB. So that's awesome. 

What about for salmon? 

Looking back at the Salmon_quant.sh script error messages, I noticed this line: "/projects/lsf_spool/1690315742.28019259.shell: line 13: 15099 Segmentation fault" as the final error that halted the script I think. 

This was for job 28019259:
```{bash}
#!/bin/bash
#BSUB -J Pcli_transcriptome_salmon_index
#BSUB -q bigmem
#BSUB -P and_transcriptomics
#BSUB -n 16
#BSUB -o Pcli_transcriptome_salmon_index%J.out
#BSUB -e Pcli_transcriptome_salmon_index%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

${and}/programs/salmon-1.5.2_linux_x86_64/bin/salmon index -t ${and}/genomes/Pcli/clean_Pcli_transcriptome_final.fasta -i ${and}/genomes/Pcli/Pcli_transcriptome_index
```

Exited with exit code 139.

Resource usage summary:

    CPU time :                                   252.31 sec.
    Max Memory :                                 1551 MB
    Average Memory :                             503.75 MB
    Total Requested Memory :                     24000.00 MB
    Delta Memory :                               22449.00 MB
    Max Processes :                              4
    Max Threads :                                7


In looking back at the salmon_index codes, I think the reason why salmon quant might not be working is because salmon_index didn't actually finish, it exited because of a fatal error. 

Or wait I'm getting confused because I submitted two jobs with the exact same indexing script above, but one to bigmem and one to general. Based on the emails, job 28019257 (the general/-n 8 node job) successfully completed, but job 28019259 (bigmem -n16) exited with a fatal error. And because the directory paths for the output files of the index were the exact same, and job 28019259 was started after job 28019257, the bigmem code overwrote the general one.

So let's try resubmitting the general one and just letting it run for awhile on its own.

I think the job holds for awhile in the queue because it waits for enough available RAM/memory before running? Idk.

I'm also going to run the job in the debug queue so I can see what it looks like.

So running it in the debug queue, it completely ran and it looks like it ran successfully and finished building the index. So I'll just wait for the general job I submitted to complete.

It worked!

I am now running salmon_quant.sh and I also submitted it to the debug queue and it seems to be working:
```{bash}
#!/bin/bash
#BSUB -J Pcli_transcriptome_salmon_quant
#BSUB -q bigmem
#BSUB -P and_transcriptomics
#BSUB -n 16
#BSUB -o Pcli_transcriptome_salmon_quant%J.out
#BSUB -e Pcli_transcriptome_salmon_quant%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/AS_pipeline/3_trimmed_fastq_files/Pcli_fastq_files"

data=($(ls *.gz))

for sample in ${data[@]} ;

do \
${and}/programs/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${and}/genomes/Pcli/Pcli_transcriptome_index -l U -r ${sample} --validateMappings -o ${and}/Ch2_temperaturevariability2023/AS_pipeline/4_Pcli_specific/salmon_quant_files ; \
done
```

Ok, so the code worked and I got the resulting file ("quant.sf"), however when I read the [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#transcript-abundance-files-and-tximport-tximeta) for how to import this into R, I see this:

<img width="926" alt="Screen Shot 2023-07-28 at 10 06 59 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/1727efa1-77ae-45c8-9d0b-bc24cbcdb71f">

When I look at the multiQC report of the trimmed fasta files, I see that all the Pcli samples "failed" the GC content thing. So I'm not sure if this means there is a bias or not, but in the [Salmon tutorial for adding -gcbias](https://salmon.readthedocs.io/en/latest/salmon.html#gcbias) it says that adding this flag when doing the quantification step "does not impair quantification for samples without GC bias, it just takes a few more minutes per sample". " For samples with moderate to high GC bias, correction for this bias at the fragment level has been shown to reduce isoform quantification errors."

So I think I'll rerun the quant code with that flag just to be safe.


<img width="1100" alt="Screen Shot 2023-07-28 at 10 07 37 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/e0d0e287-b03b-4bf1-880e-2360ed307c93">


```{bash}
#!/bin/bash
#BSUB -J Pcli_transcriptome_salmon_quant
#BSUB -q bigmem
#BSUB -P and_transcriptomics
#BSUB -n 16
#BSUB -o Pcli_transcriptome_salmon_quant%J.out
#BSUB -e Pcli_transcriptome_salmon_quant%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/AS_pipeline/3_trimmed_fastq_files/Pcli_fastq_files"

data=($(ls *.gz))

for sample in ${data[@]} ;

do \
${and}/programs/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${and}/genomes/Pcli/Pcli_transcriptome_index -l U -r ${sample} --validateMappings --gcBias --reduceGCMemory -o ${and}/Ch2_temperaturevariability2023/AS_pipeline/4_Pcli_specific/salmon_quant_files ; \
done
```

