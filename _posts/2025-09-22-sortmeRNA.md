---
layout: post
title: Using SortmeRNA
date: '2025-09-22'
categories: Analysis, Processing
tags: [reciprocal transplant, PSTR, sortmerna, Pegasus, ch3_reciprocaltransplat]
---

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


### Sept 22 update

Just checked the log files from the most recent sortmerna run, and now there's some fun new errors.

```{bash}
[write_a_read:215] ESC[0;31mERRORESC[0m: Failed deflating readstring:
[write_a_read:215] ESC[0;31mERRORESC[0m: zlib status: -1
[append:328] ESC[0;31mERRORESC[0m: Failed deflating readstring:
[append:328] ESC[0;31mERRORESC[0m: zlib status: -1
```

ChatGPT gave some bs reasons why but I honestly am just starting to give up on this program. 

Ok, I discovered that a subset of samples did work (24 out of 96):

<img width="878" height="408" alt="Screenshot 2025-09-22 at 10 03 43 AM" src="https://github.com/user-attachments/assets/5dddd384-eeed-4dd2-9067-c77cba28966b" />

Ok, I'm not sure why the rest didn't work, I'm just going to rerun them now:

```{bash}
#!/usr/bin/env bash

# Define directories
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="${and}/reciprocaltransplant"
readsdir="${projdir}/raw_seq_files/trimmed_cutadapt"
logdir="${projdir}/logs/sortmeRNA"
sortdir="${projdir}/raw_seq_files/sortmerna"
scriptdir="${projdir}/scripts/sortmerna"

# Make sure scriptdir and logdir exist
mkdir -p "${scriptdir}" "${logdir}"

# Loop through R1 trimmed reads
cd "${readsdir}" || exit

for r1 in *_R1.trimmed.fastq.gz; do
    base="${r1%_R1.trimmed.fastq.gz}"

    # Check if the SortMeRNA out folder exists and has files
    out_folder="${sortdir}/${base}/out"
    if [ -d "$out_folder" ] && [ "$(ls -A "$out_folder" 2>/dev/null)" ]; then
        echo "Skipping ${base} (output already exists)"
        continue
    fi

    # Write the LSF job script
    cat <<EOF > "${scriptdir}/${base}_sortmeRNA.job"
#!/usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${base}_sortmeRNA
#BSUB -e ${logdir}/${base}_sortmeRNA.err
#BSUB -o ${logdir}/${base}_sortmeRNA.out
#BSUB -q normal
#BSUB -n 8
#BSUB -W 48:00
#BSUB -R "rusage[mem=6000]"


# Clean up old SortMeRNA outputs for this sample
rm -rf ${projdir}/raw_seq_files/sortmerna/${base}/kvdb
cd "${readsdir}"

${and}/programs/sortmerna-4.3.6-Linux/bin/sortmerna \
  --ref ${and}/programs/sortmerna-4.3.6-Linux/database/smr_v4.3_default_db.fasta \
  --reads ${readsdir}/${base}_R1.trimmed.fastq.gz \
  --reads ${readsdir}/${base}_R2.trimmed.fastq.gz \
  --aligned ${sortdir}/${base}/${base}_rRNA \
  --other ${sortdir}/${base}/${base}_non_rRNA \
  --fastx --blast --out2 --sout \
  --workdir ${sortdir}/${base}
EOF

    # Submit the job
    bsub < "${scriptdir}/${base}_sortmeRNA.job"
done

```

Ok, this didn't work again. This time the error was "std::out_of_range / stoi". Apparently, ChatGPT is saying that version 4.3.6 of SortMeRNA is known to have a bug to cause crashes if there are any differences between reads 1 and 2 that it can't parse, I guess even if the length differ by 1 base pair. I checked my MultiQC report and I don't think that's the issue. But ChatGPT recommended removing these flags `--blast, --out2, --sout` and trying again. Otherwise, I may need to downgrade to a previous version of SortMeRNA. 

## Sept 22

Ok, I think now what's happening is I maxed out the space quota on my Pegasus account. I was at 2 T of files. I deleted the Trim Galore folder because I am using the cutadapt trimmed version files now, and that freed up ~500 GB of space. Now I'm re-running sortmeRNA, but I'm having it also delete all the subfolders for the samples that didn't work, because there could potentially be issues with having partial files written already in there. 

## Sept 23

I discovered this [open thread on the sortmerna GitHub](https://github.com/sortmerna/sortmerna/issues/379) that has the error I keep getting. It looks like there is a bug. If the file is larger than 2 GB, it stops working. 

Ok, I discovered anything greater than 20 million reads is not working in sortmerna. This equates to around 2 GB (see screenshot). Also, why is R2 larger than R1 if they have the same amount of reads?

<img width="1000" height="444" alt="Screenshot 2025-09-23 at 10 31 10 AM" src="https://github.com/user-attachments/assets/510a5f70-a14e-4239-8528-a01bac673133" />

Here are the other threads in sortmerna issues on GitHub:
- https://github.com/sortmerna/sortmerna/issues/429
- https://github.com/sortmerna/sortmerna/issues/445
- https://github.com/sortmerna/sortmerna/issues/421
- https://github.com/sortmerna/sortmerna/issues/419
- https://github.com/sortmerna/sortmerna/issues/424
- https://github.com/sortmerna/sortmerna/issues/340
- https://github.com/sortmerna/sortmerna/issues/305
- https://github.com/sortmerna/sortmerna/issues/405
- https://github.com/sortmerna/sortmerna/issues/326
  
The last few are more related to resource allocation on HPCs.

So far, the only suggestion I see is to update to use sortmerna v4.3.7. 

Someone talked about running it in parallel on several split up files: https://github.com/sortmerna/sortmerna/issues/413 

On this thread they seem to have made a good script for doing the split: https://github.com/sortmerna/sortmerna/issues/336 

One flag i saw that could be helpful is the "--no-best" function, since I'm just trying to see if any sequence aligns to any rRNA, not what specific rRNA it is (see this thread: https://github.com/sortmerna/sortmerna/issues/437).

Also I saw this flag get suggested: -threads 8

And, I'm going to put this job back in bigmem instead of general because that could be an issue too. 




