---
layout: post
title: redo Fastp for Ch2 samples
date: '2023-08-08'
categories: coding
tags: [coding, ch2_tempvariability]
---

I need to redo the Fastp step because in comparing the fastqc reports of the raw reads versus the "trimmed" reads, they are the exact same and it appears that nothing was trimmed. I also think Fastp didn't work correctly because there is only one fastp.html and fastp.json file for the entire list of samples, and it corresponds to Pcli-148 (the last sample of the batch). So I think it overwrote things and maybe didn't process correctly.

I submitted this to Pegasus:

```{bash}
#!/bin/bash
#BSUB -J trim_fastp.sh
#BSUB -q bigmem
#BSUB -n 8
#BSUB -R "rusage[mem=10000]"
#BSUB -P and_transcriptomics
#BSUB -o trim_fastp_%J.out
#BSUB -e trim_fastp_%J.err
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

module load fastqc/0.10.1

for sample in ${and}/Ch2_temperaturevariability2023/fastq_rawreads/*.fastq.gz ;
do \
${and}/programs/fastp --in1 ${sample} --out1 ${and}/Ch2_temperaturevariability2023/fastp_processed/clean.${sample} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50
fastqc ${and}/Ch2_temperaturevariability2023/fastp_processed/clean.${sample} -o ${and}/Ch2_temperaturevariability2023/trimmed_qc/ ; \
done

echo "Read trimming of adapters complete."

multiqc ${and}/Ch2_temperaturevariability2023/trimmed_qc/

mv multiqc_report.html trimmed_qc/
mv multiqc_data trimmed_qc/

echo "Cleaned MultiQC report generated."
```

But what I don't understand in the code is how to name the fastp.html files so they don't get overwritten? 

Ok I added the flag `-h report_${sample}.html` to see if it would generate a report for each sample and not overwrite them.

```{bash}
#!/bin/bash
#BSUB -J trim_fastp.sh
#BSUB -q bigmem
#BSUB -n 8
#BSUB -R "rusage[mem=10000]"
#BSUB -P and_transcriptomics
#BSUB -o trim_fastp_%J.out
#BSUB -e trim_fastp_%J.err
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

module load fastqc/0.10.1

for sample in ${and}/Ch2_temperaturevariability2023/fastq_rawreads/*.fastq.gz ;
do \
${and}/programs/fastp --in1 ${sample} --out1 ${and}/Ch2_temperaturevariability2023/fastp_processed/clean.${sample} -h report_${sample}.html --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50
fastqc ${and}/Ch2_temperaturevariability2023/fastp_processed/clean.${sample} -o ${and}/Ch2_temperaturevariability2023/trimmed_qc/ ; \
done

echo "Read trimming of adapters complete."

multiqc ${and}/Ch2_temperaturevariability2023/trimmed_qc/

mv multiqc_report.html trimmed_qc/
mv multiqc_data trimmed_qc/

echo "Cleaned MultiQC report generated."
```

We'll see if this works.

**Update**: the job failed because I'm dumb and the folders were improperly named. trying to run this now:

```{bash}
#!/bin/bash
#BSUB -J trim_fastp.sh
#BSUB -q bigmem
#BSUB -n 8
#BSUB -R "rusage[mem=10000]"
#BSUB -P and_transcriptomics
#BSUB -o trim_fastp_%J.out
#BSUB -e trim_fastp_%J.err
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

module load fastqc/0.10.1

for sample in ${and}/Ch2_temperaturevariability2023/1_fastq_rawreads/*.fastq.gz ;
do \
${and}/programs/fastp --in1 ${sample} --out1 ${and}/Ch2_temperaturevariability2023/fastp_processed/clean.${sample} -h report_${sample}.html --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50
fastqc ${and}/Ch2_temperaturevariability2023/fastp_processed/clean.${sample} -o ${and}/Ch2_temperaturevariability2023/trimmed_qc_files/ ; \
done

echo "Read trimming of adapters complete."

multiqc ${and}/Ch2_temperaturevariability2023/trimmed_qc_files/

mv multiqc_report.html trimmed_qc_files/
mv multiqc_data trimmed_qc_files/

echo "Cleaned MultiQC report generated."
```

**08.14.23 update**: 

So I actually remade this entire script because it was really slow and also it didn't do what I wanted. It did rename the .html files to each have a sample name, but that's not actually helpful for running multiqc on the fastp files (it didn't show anything). What I really needed to do was add a flag for renaming the .json files so they wouldn't overwrite themselves. So I added that flag and also made a loop job script, so it submitted a separate job for each sample to trim in parallel.

This is what the script I submitted was:
```{bash}
#BSUB -u and128@miami.edu

and="/scratch/projects/and_transcriptomics"
cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/1_fastq_rawreads"

data=($(ls *.gz))

for sample in ${data[@]} ;
do \
echo '#!/bin/bash' > "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
echo '#BSUB -q bigmem' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
echo '#BSUB -J '"${sample}"_fastp'' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
echo '#BSUB -o '"${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"$sample"_fastp%J.out'' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
echo '#BSUB -e '"${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"$sample"_fastp%J.err'' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
echo '#BSUB -n 8' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
echo '#BSUB -N' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
echo 'module load fastqc/0.10.1' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
echo '/scratch/projects/and_transcriptomics/programs/fastp '--in1 "${sample}" --out1 "${and}"/Ch2_temperaturevariability2023/fastp_processed/clean."${sample}" -h report_"${sample}".html -j "${sample}"_fastp.json --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50'' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
echo 'fastqc '"${and}"/Ch2_temperaturevariability2023/fastp_processed/clean."${sample}" -o "${and}"/Ch2_temperaturevariability2023/trimmed_qc_files/'' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
bsub < "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job ; \
done
```

As you can see in this line, I added the -j flag for renaming .json: `echo '/scratch/projects/and_transcriptomics/programs/fastp '--in1 "${sample}" --out1 "${and}"/Ch2_temperaturevariability2023/fastp_processed/clean."${sample}" -h report_"${sample}".html -j "${sample}"_fastp.json --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50'' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job`

So, the results I got were a bunch of .json files, which is good. But now I need to move them all to their own directory so I can run multiqc on just those files and not anything else.


Ok, I moved all those into a subdirectory then ran `multiqc .` there, and I finally get a report that says how many millions of reads I have post-filtering / trimming!!!

<img width="1110" alt="Screen Shot 2023-08-14 at 10 26 42 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/880a9cf9-400c-4ae5-afc5-d29fb44e355e">

<img width="1102" alt="Screen Shot 2023-08-14 at 10 27 01 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/6db0cc65-8eee-46b4-9931-f231c90ca1ea">


Now the next question is, are these trimmed fastq files the same as the originally trimmed ones? They should be because I didn't actually change any of the fastp parameters, I was just generating a multiqc report of the wrong files. I'll look at the file sizes and see if they're vastly different.

**1. original trimmed files (take_1)**:

<img width="532" alt="Screen Shot 2023-08-14 at 10 43 26 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/ad940cd8-af20-437a-a6c9-d3a410267b85">

<img width="539" alt="Screen Shot 2023-08-14 at 10 44 21 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/ca30bd37-60b9-4961-ac1f-cfa9098e65aa">

**2. recent fastp trimmed files (take_2)**:

<img width="528" alt="Screen Shot 2023-08-14 at 10 45 32 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/53430ee1-a981-440b-87a2-95f239f8866e">

<img width="535" alt="Screen Shot 2023-08-14 at 10 45 50 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/70cf5a31-007f-4611-ac18-092d84998492">

File sizes are exactly the same so that's good. I also checked the Million read numbers on the first round of fastqc/multiqc reports and those match the recent reports from the fastp files. So overall nothing changed, I just didn't have the right columns in the multiqc reports to show that adapters were trimmed and bases were filtered. 
