---
layout: post
title: Notes on Pstrigosa bioinformatics data analysis
date: '2025-08-20'
categories: Analysis
tags: [reciprocal transplant, coding, Pegasus, rRNA, de novo transcriptome assembly]
---

Dr. Brad Weiler's new publication on diel transcriptomics of *Pseudodiploria strigosa* was just published on bioRxiv [here](https://www.biorxiv.org/content/10.1101/2025.08.06.668741v1), and his pipeline is available now on GitHub [here](https://github.com/delCampoLab/diel_pstr/blob/main/pstr_diel_rnaseq_submit.Rmd). 

I want to go through his steps and follow his code since he constructed a *de novo*  transcriptome for *P. strigosa*.

The sequencing he performed was "Poly-A capture/enrichment methodology by Novogene (Beijing, China) on the
158 Illumina NovaSeq PE 150bp. The methodology follows a standard workflow of fragmentation, reverse
159 transcription, cDNA synthesis, end-repair and A-tailing."

3' RNA-Seq doesn't necessarily equate to Poly-A selection/capture/enrichment. 3' RNA-Seq can be done with no selection, meaning that it would just sequence the 3'ends of RNA molecules but isn't targeting molecules with polyA tails. Poly-A selection/capture/enrichment is performed to enrich mRNA molecules over other types or RNA (typically rRNA) that dominate total RNA content. 

The question is, did my reciprocal transplant reads get any sort of selection? Looking back through my notes, I don't see anything about Poly-A selection.

## His pipeline: 
1. **FastQC**
2. Adapter trimming using **CutAdapt** wrapper script TrimGalore
3. Ribosomal RNA partitioned from total RNA using **SortMeRNA** and default fasta reference database
4. Remaining RNA mapped to index of genomes (Ostreobium quekettii,
Symbiodinium sp. 57, Breviolum sp. 58, Cladocopium sp. 57, Durusdinium sp. 59 183 , Eimeria tenella, Emitis
184 houghton, Sarcocystis neurona, Toxoplasma gondii, Malassezia restricta, Nitzschia putrida, Gracilaria
gracilis, and Uronemita sp.) to isolate non-coral RNA reads using **BBsplit**
5. Unmapped reads concaenated for *de novo* transcriptome assembly using **Trinity**
6. Decontamination using Standard PlusPF reference database in **Kraken2**
7. Reads aligned to transcriptome using **Salmon**
8. **STAR** used to align *Breviolum* binned reads
9. Gene transcripts were first predicted into longest open reading frames (ORFs) using **transDecoder**
10. Annotations using **EggNOG-mapper**
11. Gene ontology (GO) annotations supplemented using GoPredSim ProtT5 protein language model 65 195 derived topGO terms using **FANTASIA**

### Using SortMeRNA for removing rRNA
- Looking online, there seems to be some thoughts that removing rRNA reads is not required when Poly-A selection is performed during sequencing, as preferentially selecting polyA tails automatically creates a bias towards mRNA (see [this article](https://www.biostars.org/p/419845/) and [this article]()).
- However, some suggest that removing rRNA would save computation time and lead to a cleaner assembly ([see here](https://www.researchgate.net/post/Should_I_remove_rRNAs_from_transcriptome_data_while_moving_on_to_Denovo_assembly)).
- This [tutorial/paper](https://trhvidsten.com/docs/Delhomme-EpiGeneSysProtocol2015.pdf) is a good resource - it says that if GC content is abnormal, it could be a sign of high rRNA content (GC content of over 50% is typical for rRNA) and then SortMeRNA should be used. They also said "if there is any doubt, this step should be performed."

Based on the multiQC report of my trimmed reads, it looks like the average GC content is not normally distributed and is higher for some samples. To be safe, I'll follow Brad's code and run **SortMeRNA** on it.

<img width="959" height="528" alt="Screenshot 2025-08-20 at 1 54 00 PM" src="https://github.com/user-attachments/assets/82bfdc4b-c67a-42c8-bafb-a946812b8e2f" />

#### Update

So I ran SortMeRNA. Now what? I don't want to do the BBSplit step because I'm not looking to see what aligns to members of the microbiome. I think I just want to skip to Trinity. But what files from SortMeRNA do I use to do that?

When I look at Brad's next steps for BBSplit and Trinity, it looks like he renamed his fastq files from SortMeRNA. 

The results of SortMeRNA are directories for each sample. 
<img width="351" height="31" alt="Screenshot 2025-09-13 at 5 47 17 PM" src="https://github.com/user-attachments/assets/ceaf224e-c0bf-4bf2-ba30-9ab45732cd43" />

According to the [SortMeRNA user manual](https://sortmerna.readthedocs.io/en/latest/manual4.0.html#usage): 
- "kvdb" = "key-value datastore for alignment results"
- idx = "index database"
- "out" = "output files like aligned.blast"

It doesn't mention the other folder I have, which is "readb", and that's the one I found .fq.gz files in, which look like what I would work with?

Ok I found on [this website](https://hpc.nih.gov/apps/sortmerna.html) more info: it says that "readb" means "pre-processed reads/index". 

ChatGPT said that my code isn't working because I didn't assign it a name for the aligned versus unaligned reads. Idk... Working on trying a new code it suggested now:

```{bash}
#!/usr/bin/env bash

# Define project directories and paths
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="${and}/reciprocaltransplant"
readsdir="${projdir}/raw_seq_files/trimmed_cutadapt"
logdir="${projdir}/logs/sortmeRNA"

cd "${readsdir}"

# Loop through R1 files
for r1 in *_R1.trimmed.fastq.gz; do
  # Extract sample base name (e.g., "Pstr-Dec2022-202")
  base="${r1%_R1.trimmed.fastq.gz}"

  # Submit the job directly
  bsub <<EOF
#!/usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${base}_sortmeRNA
#BSUB -e ${logdir}/${base}_sortmeRNA.err
#BSUB -o ${logdir}/${base}_sortmeRNA.out
#BSUB -q bigmem
#BSUB -n 10

cd "${readsdir}"

${and}/programs/sortmerna-4.3.6-Linux/bin/sortmerna \\
  --ref ${and}/programs/sortmerna-4.3.6-Linux/database/smr_v4.3_default_db.fasta \\
  --reads ${readsdir}/${base}_R1.trimmed.fastq.gz \\
  --reads ${readsdir}/${base}_R2.trimmed.fastq.gz \\
  --aligned ${base}_rRNA \\
  --other ${base}_non_rRNA \\
  --fastx \\
  --out2 \\
  --sout \\
  --log \\
  --workdir ${projdir}/raw_seq_files/sortmerna/${base}
EOF
done
```

**Update**: this script sat in the job queue for 2 days. I don't think it will work. 

I also realized that I don't think it's using the files that sortmeRNA already made. I'm having chatGPT write a few versions to try to see what will work. 


**Sept 17 2025 Update**:
- I realized in the original error files for sortmeRNA, I got "run time limit reached" flags for each job. So I think what happened was that they all timed out at around 6 hours. In an unrelated note, but may be related, the Pegasus technicians told me that I needed to start using pegasus2 rather than pegasus. I discovered this when I tried to check the status of my jobs, and the LSF "daemon" wasn't responding.

- So I have done that, and I'm basically trying to re-run the original script that I got to work, but that timed out. I added the specific file names for --aligned and --other in the hopes that it would help. I also added a longer time limit request in the job submission, so hopefully that helps.

The job I ran most recently is called "6_sortmerna.sh" on Pegasus and looks like this:

```{bash}
#!/usr/bin/env bash

# Define project directories and paths
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="${and}/reciprocaltransplant"
readsdir="${projdir}/raw_seq_files/trimmed_cutadapt"
logdir="${projdir}/logs/sortmeRNA"

cd "${readsdir}"

# Loop through R1 files
for r1 in *_R1.trimmed.fastq.gz; do
  # Extract sample base name (e.g., "Pstr-Dec2022-202")
  base="${r1%_R1.trimmed.fastq.gz}"

# Write the job script
    cat <<EOF > "${projdir}/scripts/sortmeRNA/${base}_sortmeRNA.job"

#!/usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${base}_sortmeRNA
#BSUB -e ${logdir}/${base}_sortmeRNA.err
#BSUB -o ${logdir}/${base}_sortmeRsNA.out
#BSUB -q bigmem
#BSUB -n 10
#BSUB -W 48:00
#BSUB -R "rusage[mem=6000]"

cd "${readsdir}"

${and}/programs/sortmerna-4.3.6-Linux/bin/sortmerna --ref ${and}/programs/sortmerna-4.3.6-Linux/database/smr_v4.3_default_db.fasta \
--reads ${readsdir}/${base}_R1.trimmed.fastq.gz \
--reads ${readsdir}/${base}_R2.trimmed.fastq.gz \
--aligned --other --fastx --blast --out2 --sout \
--workdir ${projdir}/raw_seq_files/sortmerna/${base}

EOF

    # Submit the job
    bsub < "${projdir}/scripts/sortmerna/${r1}_sortmerna.job"
done
```

but that last line about submitting the jobs doesn't work, so I'll need to manually submit all the jobs.

I ran this in the command line: 

```{bash}
for f in *.job; do
    bsub < "$f"
done
```

Ok, the "bigmem" queue ones never submitted and were still pending a day later, so I used `bjobs -p` to see why, and it said this: "Not enough processors to meet the job's spanning requirement: 1 host;"

So I changed the parameters to "normal" queue and -8 processors, and it ran. But then all the jobs exited out because of the following issue: "Please, ensure the directory "/scratch/projects/and_transcriptomics/reciprocaltransplant/raw_seq_files/sortmerna/Pstr-Nov2023-198/kvdb" is Empty prior running 'sortmerna'".

According to chatGPT, having a pre-filled "kvdb" folder is the only one that sortmerna cares about and will refuse to rewrite over.

So now I'm running this script to erase those subfolders for all samples in advance, then trying again.

```{bash}
#!/usr/bin/env bash

# Define project directories and paths
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="${and}/reciprocaltransplant"
readsdir="${projdir}/raw_seq_files/trimmed_cutadapt"
logdir="${projdir}/logs/sortmeRNA"

cd "${readsdir}"

# Loop through R1 files
for r1 in *_R1.trimmed.fastq.gz; do
  # Extract sample base name (e.g., "Pstr-Dec2022-202")
  base="${r1%_R1.trimmed.fastq.gz}"

# Write the job script
    cat <<EOF > "${projdir}/scripts/sortmeRNA/${base}_sortmeRNA.job"

#!/usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${base}_sortmeRNA
#BSUB -e ${logdir}/${base}_sortmeRNA.err
#BSUB -o ${logdir}/${base}_sortmeRsNA.out
#BSUB -q normal
#BSUB -n 8
#BSUB -W 48:00
#BSUB -R "rusage[mem=6000]"


# Clean up old SortMeRNA outputs for this sample
rm -rf ${projdir}/raw_seq_files/sortmerna/${base}/kvdb

cd "${readsdir}"

${and}/programs/sortmerna-4.3.6-Linux/bin/sortmerna --ref ${and}/programs/sortmerna-4.3.6-Linux/database/smr_v4.3_default_db.fasta \
--reads ${readsdir}/${base}_R1.trimmed.fastq.gz \
--reads ${readsdir}/${base}_R2.trimmed.fastq.gz \
--aligned --other --fastx --blast --out2 --sout \
--workdir ${projdir}/raw_seq_files/sortmerna/${base}

EOF

    # Submit the job
    bsub < "${projdir}/scripts/sortmerna/${r1}_sortmerna.job"
done
```


