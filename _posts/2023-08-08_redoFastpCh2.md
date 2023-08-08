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
