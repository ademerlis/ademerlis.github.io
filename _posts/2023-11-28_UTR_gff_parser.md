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

Ok I figured out the issue with the code, it was having an issue with grep "#" and it needed to be grep '#' , then everything worked. But I'm trying to compare the Acer_parsed.gtf file to the Acer_parsed_ext_5000.gtf, and it's really hard to tell if they're different. 

First, here's the code that I got to work:

```{bash}
#!/usr/bin/env bash
#BSUB -P and_transcriptomics
#BSUB -e /scratch/projects/and_transcriptomics/genomes/Acer_2023/UTR_add_extend_GTF-main/UTR_add.err
#BSUB -o /scratch/projects/and_transcriptomics/genomes/Acer_2023/UTR_add_extend_GTF-main/UTR_add.out

module load python/3.8.7
cd /scratch/projects/and_transcriptomics/genomes/Acer_2023/UTR_add_extend_GTF-main/
GTF_FILE="Acer_parsed.gtf"

mkdir temp_files;
mkdir results;
grep '#' ${GTF_FILE} > ${GTF_FILE}_hashtag_lines
# split input gtf into fw and rev
awk -F'\t' '$7=="+" {print $0}' ${GTF_FILE} > ${GTF_FILE}_fw.gtf
awk -F'\t' '$7=="-" {print $0}' ${GTF_FILE} > ${GTF_FILE}_rev.gtf

cd /scratch/projects/and_transcriptomics/genomes/Acer_2023/UTR_add_extend_GTF-main/ ;
module load python/3.8.7 ;
python3 Extend_by_threshold_fw.py --input ${GTF_FILE}_fw.gtf --threshold 5000 > temp_files/${GTF_FILE}_fw.gtf ;
python3 Extend_by_threshold_rev.py --input ${GTF_FILE}_rev.gtf --threshold 5000 > temp_files/${GTF_FILE}_rev.gtf ;
cat temp_files/${GTF_FILE}_fw.gtf > temp_files/unsorted ;
cat temp_files/${GTF_FILE}_rev.gtf >> temp_files/unsorted ;
sort -k1,1V -k4,4n -k5,5rn temp_files/unsorted > temp_files/sorted ;
cat ${GTF_FILE}_hashtag_lines > results/${GTF_FILE}_ext_5000.gtf ;
cat temp_files/sorted >> results/${GTF_FILE}_ext_5000.gtf
```

Then when i look at the results file, it's the exact same size as the Acer_parsed.gtf file.

But, when i run `diff -y ../Acer_parsed.gtf Acer_parsed.gtf_ext_5000.gtf , I get a bunch of lines that are different. It looks like in places there were gaps, it created exons, genes, and transcripts that fill in the gaps. Which is confusing because shouldn't it be labeling them as UTRs? 

**November 29 update**:

Ok I think I finally understand what this code is doing thanks to Nick explaining it! 

Essentially the gtf_advanced_parser.py is the script that changes the annotations in column 3 of the gtf file.

Here's an example for one gene of what the original file looks like (GCA_032359415.1_NEU_Acer_K2_genomic.gtf):

![Screen Shot 2023-11-29 at 1 57 10 PM](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/a582039b-34ba-435a-8b36-3671d2ef9cb6)

Versus the file after running gtf_advanced_parser.py:

![Screen Shot 2023-11-29 at 1 57 43 PM](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/f0ec8e42-67b1-4085-bbf7-9716576bc4df)


And then the same gene after the process_them_threshold5000.sh job is run: 

![Screen Shot 2023-11-29 at 1 58 09 PM](https://github.com/ademerlis/ademerlis.github.io/assets/56000927/3490a8db-cb74-4d72-8a97-6ad3bde91d78)


It looks like what this code overall is doing is removing "CDS" and start/stop codons, just keeping exons, genes, and transcripts, and extending the lengths of the gene, transcript and exon so that they match up to the next adjacent exon. Essentially removing gaps in the gtf file. It's artifically making the genes longer. So I guess this hypothetically accounts for UTRs by replacing any gaps with just longer exons of the same gene? But the thing I was getting confused by was that nothing was strictly labeled as a UTR in the gff or gtf file. 

Ok so now that I'm sure that this UTR gtf parser thing worked (in terms of adding numbers to extend the exons/genes/transcripts), I can move on to working on the alignment step with the new Acer genome. 
