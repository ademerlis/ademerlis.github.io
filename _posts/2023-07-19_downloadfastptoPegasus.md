---
layout: post
title: download fastp to Pegasus
date: '2023-07-19'
categories: coding
tags: [coding, temperaturevariability2023]
---

The next step of the pipeline is to run fastp to conduct trimming and cleaning of sequence files.

[Zoe runs this code](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/ZD_Heron-Pdam-gene-expression.md):

```{bash}
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --error=../"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=../"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/zdellaert/Pdam-TagSeq/raw_data

# load modules needed
module load fastp/0.19.7-foss-2018b
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.9-intel-2020a-Python-3.8.2

#make qc output folder
mkdir /data/putnamlab/zdellaert/Pdam-TagSeq/trimmed_qc/

#make processed folder for trimmed reads
mkdir /data/putnamlab/zdellaert/Pdam-TagSeq/processed/

# Make an array of sequences to trim
array1=($(ls *.fastq.gz)) 

# fastp loop; trim the Read 1 TruSeq adapter sequence; trim poly x default 10 (to trim polyA) 

for i in ${array1[@]}; do
 fastp --in1 ${i} --out1 /data/putnamlab/zdellaert/Pdam-TagSeq/processed/clean.${i} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50 
        fastqc /data/putnamlab/zdellaert/Pdam-TagSeq/processed/clean.${i} -o /data/putnamlab/zdellaert/Pdam-TagSeq/trimmed_qc/ # fastqc the cleaned reads
done 

echo "Read trimming of adapters complete." $(date)

# Quality Assessment of Trimmed Reads

multiqc /data/putnamlab/zdellaert/Pdam-TagSeq/trimmed_qc/ #Compile MultiQC report from FastQC files 

mv multiqc_report.html trimmed_qc/ #move output files to the QC directory
mv multiqc_data trimmed_qc/ #move output files to the QC directory

echo "Cleaned MultiQC report generated." $(date)
```

I need to download a local/scratch space source of the fastp program because Pegasus doesn't have that module installed. I've tried a couple things so far and neither worked:
1. cloning source: "git clone https://github.com/OpenGene/fastp.git" (then you need to make and install and then also install all the dependencies manually - I think this will work if I do it all, it just takes awhile and idk how many dependencies I need to install so I want it to do it all for me)
2. module load anaconda on pegasus then do conda install: "conda install -c bioconda fastp" (it then asks me to update conda and i can't do that for pegasus because it's a shared module and i need admin access)

[methods for installing fastp](https://github.com/OpenGene/fastp)

I ran this in pegasus in my programs folder and i think it worked:
```{bash}
# download the latest build
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp
```

the only way to test it is to try running a script with it now.

```{bash}
#!/bin/bash
#BSUB -J trim_qc.sh
#BSUB -n 8
#BSUB -R "rusage[mem=10000]"
#BSUB -P and_transcriptomics
#BSUB -o trim_qc_%J.out
#BSUB -e trim_qc_%J.err
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

# load modules needed
module load fastqc/0.10.1

#make qc output folder
mkdir ${and}/Ch2_temperaturevariability2023/trimmed_qc/

#make processed folder for trimmed reads
mkdir ${and}/Ch2_temperaturevariability2023/processed/

# Make an array of sequences to trim
array1=($(ls *.fastq.gz)) 

# fastp loop; trim the Read 1 TruSeq adapter sequence; trim poly x default 10 (to trim polyA) 

cd ${and}/Ch2_temperaturevariability2023/fastq_rawreads/

for i in ${array1[@]}; do
 ${and}/programs/fastp --in1 ${i} --out1 ${and}/Ch2_temperaturevariability2023/processed/clean.${i} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50 
        fastqc ${and}/Ch2_temperaturevariability2023/processed/clean.${i} -o ${and}/Ch2_temperaturevariability2023/trimmed_qc/ # fastqc the cleaned reads
done 

echo "Read trimming of adapters complete." $(date)

# Quality Assessment of Trimmed Reads

multiqc ${and}/Ch2_temperaturevariability2023/trimmed_qc/ #Compile MultiQC report from FastQC files 

mv multiqc_report.html trimmed_qc/ #move output files to the QC directory
mv multiqc_data trimmed_qc/ #move output files to the QC directory

echo "Cleaned MultiQC report generated." $(date)
```

Let's give it a try. 

While it's running: note that in all their codes, they specify this adapter sequence: --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

That's because it corresponds to the [Illumina TruSeq single index and CD index based kits read 1](https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314).

