---
layout: post
title: Tximport for salmon quant Pcli
date: '2023-08-03'
categories: coding
tags: [coding, ch2_tempvariability]
---

So I had a lot of issues in trying to import the "quant.sf" files into R so that I could generate a counts matrix of transcripts. 

First, each quant.sf file for each sample is massive (~40 MB) so trying to upload 48 of those files to GitHub (which my GitHub repo is locally stored on my iCloud drive which only has 200 GB of space) took too long. So I stored all the salmon_quant_files from Pegasus onto my OneDrive instead. 

Then, I had issues following this [tutorial](http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#Salmon). But this is what I ended up getting to work after trial and error AND asking ChatGPT every step of the way what to do:

```{r}
library("tximport") #requires BiocManager
library("readr")
library("tximportData")
library("readxl")
library(tidyverse)
library(DESeq2)

main_dir <- file.path("~/OneDrive - University of Miami/NOAA ERL/stress hardening 2022/gene expression/Pclivosa/salmon_quant_files/")

#make a list of subdirectories
sub_dir <- list.dirs(main_dir, recursive = F)

#initialize empty list to store all "quant.sf" file paths
quant_sf_files <- list()

#Loop through each subdirectory and get file path for "quant.sf"
for (subdir in sub_dir) {
  sample_name <- tools::file_path_sans_ext(basename(subdir))
  quant_sf_path <- file.path(subdir, "quant.sf")
  quant_sf_files[[subdir]] <- normalizePath(quant_sf_path)
}

# Extract the file paths from the named list
file_paths <- unlist(quant_sf_files)

# Check if the files exist
all_exist <- all(file.exists(file_paths))

#import quant.sf files
txi_data <- tximport(files = file_paths, type = "salmon", txOut=TRUE)
#NOTE: you need the txOut=TRUE flag if you want to summarize to transcript-level. The default for tximport is to summarize to gene-level, but we can't do that with this dataset because we are using a de novo transcriptome built via Trinity.

names(txi_data)
#abundance #counts #length #countsFromAbundance
```
Tidying counts matrix so that I can have clean sample names
```{r}
head(txi_data$counts)

colnames(txi_data$counts)

countsmatrix <- as.data.frame(txi_data$counts)

sample_names <- sub(".*/clean\\.([^_]*)_.*", "\\1", x = quant_sf_files)
sample_names <- gsub("-", "_", sample_names)

colnames(countsmatrix) <- sample_names

#round each value to the nearest integer (DESeq2 can't take decimals)
countsmatrix <- apply(countsmatrix, 2, round)

countsmatrix_df <- as.data.frame(countsmatrix)

write_csv(countsmatrix_df, "../../../results/Pcli_transcript_counts_matrix.csv")
```
making sample metadata so that I can match treatment and colony to each sample
```{r}
sample_names

sample_metadata <- as.data.frame(sample_names)

samples <- read_csv("Pcli_samples.csv")

samples %>% 
  select(colony, treatment, `Fastq file name`) %>% 
  drop_na() %>% 
  mutate(sample_names = sub("clean\\.([^_]*)_.*", "\\1", `Fastq file name`)) %>% 
  mutate(sample_names = gsub("-", "_", sample_names)) %>% 
  select(!`Fastq file name`) -> samples_tidy

column_to_rownames(samples_tidy, var = "sample_names") -> samples_tidy
```

Now I can finally make the DESeq2 object!!
```{r}
dds <- DESeqDataSetFromMatrix(countsmatrix, samples_tidy, ~ treatment)
```

I think I want to try to run 2 models though, because I should also see if Colony (genotype but we didn't genotype the corals) plays a role in gene expression (it likely does). 

Ok wait but before I create the dds object, I need to make sure all the samples had good enough alignment and high enough reads. Because samples with really low alignment rates will mess with the DEG analysis. I'm not sure if the read count (i.e. at least 4 million reads) is a necessary flag for Salmon since it's pseudoalignment. But at least for obtaining mapping rates (reads from my samples that mapped to the reference transcriptome), that needs to be >50%. 

To obtain this, I have to extract the mapping rates from the salmon_quant.log file from each sample.

Check mapping rates for each sample
```{r}
list.dirs(sub_dir)

#make a list of directories for logs
log_dir <- list.dirs(sub_dir, recursive = F)

log_files <- list.files(log_dir, pattern = "\\.log$", full.names = T)

# Initialize an empty vector to store the mapping rates
mapping_rate_lines <- character()
file_names <- character()

for (log_file in log_files) {
  # Read all lines from the log file
  log_lines <- readLines(log_file)
  
  # Find the line containing the mapping rate using the stringr package
  mapping_rate_line <- log_lines[str_detect(log_lines, "Mapping rate")]
  
  # If the line is found, extract the mapping rate value
  if (length(mapping_rate_line) > 0) {
    mapping_rate_line <- gsub("\\[.*?\\]", "", mapping_rate_line)

    # Modify the log_file name using sub
    modified_log_file <- sub(".*/clean\\.([^_]*)_.*", "\\1", log_file)
    
    file_names <-c(file_names, modified_log_file) 
    mapping_rate_lines <- c(mapping_rate_lines, mapping_rate_line)
  }
}

# Print the mapping rates
print(mapping_rate_lines)

#create data frame with file names and mapping rate lines as separate columns
mapping_rate_data <- data.frame(File_Name = file_names, Mapping_Rate_Line = mapping_rate_lines)

# Write the data frame to a tab-separated file
write.table(mapping_rate_data, file = "mapping_rates.txt", sep = "\t", row.names = FALSE)
```

All of my samples have > 50% alignment so yay. 



