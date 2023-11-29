---
layout: post
title: UTR gff parser
date: '2023-11-28'
categories: Coding
tags: [Coding]
---

I'm going to try to run [this](https://github.com/danilotat/UTR_add_extend_GTF/tree/main) again because Nick said he had no issues doing it. I think maybe if I try to run gtf_advanced_parser.py first on the genomic.gff file downloaded from NCBI, then maybe I can see if that works and I can find the "three_prime_UTR" annotations in the gff file.

When I run this, I think this works:

```{bash}
(base) [and128@login4 UTR_add_extend_GTF-main]$ python3 gtf_advanced_parser.py --input /scratch/projects/and_transcriptomics/genomes/Acer_2023/Vollmer_2023/GCA_032359415.1_NEU_Acer_K2_genomic.gtf --output Acer_parsed.gtf
We finished here. Here's some info:
Processed genes: 33794
Processed transcripts: 28059
Added 3'UTR onto forward strand: 13963
Added 3'UTR onto reverse strand: 14096
Report any error to @danilotat [ github.com/danilotat ]

Enjoy!
```

But when I run part 2 of the code, I get files with nothing in them.

```{bash}
#!/usr/bin/env bash
#BSUB -P and_transcriptomics
#BSUB -e /scratch/projects/and_transcriptomics/genomes/Acer_2023/UTR_add_extend_GTF-main/UTR_add.err
#BSUB -o /scratch/projects/and_transcriptomics/genomes/Acer_2023/UTR_add_extend_GTF-main/UTR_add.out

module load python/3.8.7
cd /scratch/projects/and_transcriptomics/genomes/Acer_2023/UTR_add_extend_GTF-main/
GTF_FILE="Acer_parsed.gtf"

# grep hashtag lines

mkdir temp_files;
mkdir results;
grep "#" ${GTF_FILE} > ${GTF_FILE}_hashtag_lines
# split input gtf into fw and rev
awk -F'\t' '$7=="+" {print $0}' ${GTF_FILE} > ${GTF_FILE}_fw.gtf
awk -F'\t' '$7=="-" {print $0}' ${GTF_FILE} > ${GTF_FILE}_rev.gtf

for threshold in "5000"; do
        cd /scratch/projects/and_transcriptomics/genomes/Acer_2023/UTR_add_extend_GTF-main ;
        module load python/3.8.7 ;
        python3 Extend_by_threshold_fw.py --input ${GTF_FILE}_fw.gtf --threshold $threshold > temp_files/${GTF_FILE}_fw.gtf ;
        python3 Extend_by_threshold_rev.py --input ${GTF_FILE}_rev.gtf --threshold $threshold > temp_files/${GTF_FILE}_rev.gtf ;
        cat temp_files/${GTF_FILE}_fw.gtf > temp_files/unsorted ;
        cat temp_files/${GTF_FILE}_rev.gtf >> temp_files/unsorted ;
        sort -k1,1V -k4,4n -k5,5rn temp_files/unsorted > temp_files/sorted ;
        cat ${GTF_FILE}_hashtag_lines > results/${GTF_FILE}_ext_by_${threshold}.gtf ;
        cat temp_files/sorted >> results/${GTF_FILE}_ext_by_${threshold}.gtf ;
        rm -rf temp_files/* ;
done

rm ${GTF_FILE}_hashtag_lines ;
rm ${GTF_FILE}_fw.gtf ;
rm ${GTF_FILE}_rev.gtf ;
```

Something is wrong in this code or it's not finding what it is supposed to, because the _fw.gtf and _rev.gtf and _hashtag_lines files are all 0 kb in size.
