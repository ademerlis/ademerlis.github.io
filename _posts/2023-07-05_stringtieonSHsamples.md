---
layout: post
title: stringtie on SH samples
date: '2023-07-05'
categories: coding
tags: [coding, temperaturevariability2023]
---

I am now going to try the stringtie scripts on the stress-hardening Acer samples and see if I get better results than the Acer CCC ones. 

I still get the same error: "/projects/lsf_spool/1688568095.27980042.shell: line 21:  : command not found"

this is the script I ran:
```{bash}
#!/bin/bash
#BSUB -J stringtie
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o stringtie%J.out
#BSUB -e stringtie%J.err
#BSUB -u and128@miami.edu
#BSUB -N

module load python/3.8.7

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Allyson_stresshardening_RNAseq/aligned"

data=($(ls *Aligned.sortedByCoord.out.bam))

for i in ${data[@]} ;

do \
/scratch/projects/and_transcriptomics/programs/stringtie-2.2.1/stringtie -G /scratch/projects/and_transc$
done
```

The next step is to merge the stringTie-produced .gtf files

```{bash}
mkdir stringtie_gtf_files
mv *.gtf stringtie_gtf_files/
ls *gtf > acerv_mergelist.txt
cat acerv_mergelist.txt
```

Now I run the merge script:

```{bash}
#!/bin/bash
#BSUB -J stringtie_mergegtf
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o stringtie_mergegtf%J.out
#BSUB -e stringtie_mergegtf%J.err
#BSUB -u and128@miami.edu
#BSUB -N

module load python/3.8.7

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Allyson_stresshardening_RNAseq/aligned/stringtie_gtf_files"

${and}/programs/stringtie-2.2.1/stringtie --merge -p 8 -G ${and}/genomes/Acer/Acerv_assembly_v1.0.gff3 -o stringtie_acerv_merged.gtf acerv_mergelist.txt
```

Next, assess assembly quality using gffcompare:
```{bash}
cd /scratch/projects/and_transcriptomics/Allyson_stresshardening_RNAseq/aligned/stringtie_gtf_files

/scratch/projects/and_transcriptomics/programs/gffcompare-0.12.6.Linux_x86_64/gffcompare -r /scratch/projects/and_transcriptomics/genomes/Acer/Acerv_assembly_v1.0.gff3 -o Acerv.merged stringtie_acerv_merged.gtf
48478 reference transcripts loaded.
  101293 query transfrags loaded
```

Now look at the summary file (Acerv.merged): 

<img width="799" alt="Screen Shot 2023-07-05 at 11 39 04 AM" src="https://github.com/ademerlis/ademerlis.github.io/assets/56000927/00aa7b4b-c412-4dd2-bcf3-47712bf62ff9">

Next, re-estimate assembly:

```{bash}
#!/bin/bash
#BSUB -J stringtie_reestimateassembly
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o stringtie_reestimateassembly%J.out
#BSUB -e stringtie_reestimateassembly%J.err
#BSUB -u and128@miami.edu
#BSUB -N

module load python/3.8.7

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Allyson_stresshardening_RNAseq/aligned"

data=($(ls *Aligned.sortedByCoord.out.bam))

for i in ${data[@]} ;

do \
${and}/programs/stringtie-2.2.1/stringtie -e -G /${and}/genomes/Acer/Acerv_assembly_v1.0.gff3 -o ${i}.merge.gtf ${i} ; \ 
done
```

Notes on re-estimating assembly (the above code doesn't make sense to me):

So I am getting confused with this step of the analysis because it doesn't make sense to me to re-estimate the assembly of the .sorted.bam samples with the same .gff3 genome file that I used in the first step of the stringtie pipeline -- i'm basically just re-doing step 1 without using the results of gffcompare or the merged gtf files. So I asked ChatGPT why it is important to re-estimate assembly after using stringtie and gffcompare, and this is what it said:

"When using StringTie and GffCompare for RNA-seq data analysis, it is common to re-estimate the assembly after running GffCompare to improve the accuracy and completeness of the transcript annotation. Here's why:

1. Initial assembly limitations: The initial assembly performed by StringTie generates a transcriptome annotation based on the input reads and reference annotation (if available). However, the initial assembly may have limitations such as incomplete or fragmented transcripts, inaccurate exon boundaries, or missed alternative splicing events.

2. GffCompare analysis: GffCompare compares the assembled transcripts with a reference annotation to identify novel isoforms, validate existing ones, and assess their compatibility with the reference. It assigns various status labels (e.g., "u" for novel isoforms, "j" for compatible junctions) to each assembled transcript based on its relationship with the reference annotation.

3. Refinement of assembly: After GffCompare analysis, re-estimating the assembly with StringTie allows the incorporation of the GffCompare results to improve the assembly quality. This step takes into account the information from GffCompare, including the novel isoforms and the compatibility of junctions, to refine the transcript structures and update their expression estimates.

4. Enhanced annotation accuracy: By re-estimating the assembly, StringTie can revise the transcript structures, adjust exon boundaries, and incorporate novel isoforms discovered by GffCompare. This process results in an enhanced and more accurate annotation of the transcriptome, which is crucial for downstream analyses such as differential expression analysis or functional interpretation.

Re-estimating the assembly after running GffCompare helps to iteratively improve the transcriptome annotation by integrating information from both the initial assembly and the reference annotation. It is an important step to obtain a comprehensive and reliable transcriptome representation for further analyses." 

This makes sense to me. You use the results of the gffcompare to re-estimate. But when I look at [Jill's code](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/RNASeq_pipeline_FL.md), I get confused because in her section for re-estimating assembly, she runs "Acerv_assembly_v1.0.gff3" as the genome for stringtie against the original .bam files, not the updated gtf file from gffcompare. So I asked ChatGPT whether if you run gffcompare and then stringtie sequentially, it "knows" that you have an updated file, and no it doesn't: 

"Apologies for the confusion in my previous response. It seems there was a misunderstanding regarding the specific workflow involving StringTie and GffCompare. Let me clarify the correct process:

1. Run StringTie: Initially, you run StringTie to assemble and annotate the transcripts based on your RNA-seq data. This step generates a GTF (Gene Transfer Format) file containing the assembled transcripts.

2. Run GffCompare: After running StringTie, you would then run GffCompare on the assembled transcripts (GTF file) from Step 1 and a reference annotation file (GTF/GFF format). GffCompare compares the assembled transcripts to the reference annotation, identifies novel isoforms, assesses their compatibility with the reference, and assigns status labels accordingly.

3. Generate updated GTF file: GffCompare generates an updated GTF file that incorporates the comparison results and assigns unique transcript IDs and class codes to each transcript. This updated GTF file contains the improved transcript structures and additional information based on the comparison with the reference annotation.

4. Use the updated GTF file: Once you have the updated GTF file generated by GffCompare, you can use it in subsequent analyses. The updated GTF file represents the refined and improved annotation of the transcriptome based on the information from GffCompare. You can use this updated GTF file for downstream analyses such as differential expression analysis or functional interpretation.

It's important to note that you don't need to re-run StringTie after GffCompare. The GffCompare step is performed on the output of StringTie to enhance the transcript annotation. The generated updated GTF file from GffCompare is used directly for further analysis.

I apologize for any confusion caused by the incorrect information provided earlier. Please let me know if you have any further questions!"

So I think I need to edit the above script for re-estimating assembly. 

```{bash}
#!/bin/bash
#BSUB -J stringtie_reestimateassembly
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o stringtie_reestimateassembly%J.out
#BSUB -e stringtie_reestimateassembly%J.err
#BSUB -u and128@miami.edu
#BSUB -N

module load python/3.8.7

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Allyson_stresshardening_RNAseq/aligned"

data=($(ls *Aligned.sortedByCoord.out.bam))

for i in ${data[@]} ;

do \
${and}/programs/stringtie-2.2.1/stringtie -e -G /${and}/Allyson_stresshardening_RNAseq/aligned/stringtie_gtf_files/stringtie_acerv_merged.gtf -o ${i}_reestimated.merge.gtf ${i} ; \ 
done
```

I think then we use the .merged_gtf files produced from this stringtie analysis for each sample to extract read count information.

(ran this in command line no job)
```{bash}
F="/scratch/projects/and_transcriptomics/Allyson_stresshardening_RNAseq/aligned/stringtie_reestimatedmerged_gtf_files/"

data=($(ls *merge.gtf))

for i in ${data[@]} ;

do \
echo "${i} $F${i}" >> sample_list_acerv.txt ; \
done

/scratch/projects/and_transcriptomics/programs/stringtie-2.2.1/prepDE.py -g gene_count_acerv_matrix.csv -i sample_list_acerv.txt
```

Then move gene counts onto local folder and you're all set woooooooooo

