---
layout: post
title: renaming trimmed files
date: '2023-07-23'
categories: coding
tags: [coding]
---

Following the fastp trimming code, all the file names ended in ".clean.processed" and that prevented them from being recognized as fastq.gz files for the next step in the pipeline. So I needed to write a for-loop to remove those suffices from all file names. 

what ended up working was:

```{bash}
files=($(ls /scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/AS_pipeline/trimmed))

for file in "${files[@]}"; do     
filename="${file%.*}";    
echo "Original: $file  |  Without Suffix: $filename";          
mv "$file" "$filename"; 
done
```

I ran this code twice (and had to redefine the files variable in between) so that it would first remove ".processed" then remove ".clean". The "%.*" part is what tells it to remove anything after the last period, including the period. 

Thanks ChatGPT!
