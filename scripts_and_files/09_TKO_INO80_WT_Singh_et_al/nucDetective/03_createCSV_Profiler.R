#!/usr/bin/env Rscript 

# Get working directory
pathSTD<-getwd()

# Load libraries
library(stringr)

# List fastq files
files<-list.files("data/raw_data/")

# Define forward/reverse reads
fwd<-grep("1.fastq.gz",files, value = T)
rev<-grep("2.fastq.gz",files, value = T)

# Define sample name (removing the file ending)
name<-fwd %>% str_sub(24,-1) %>% str_split_i("_Sac",1)

#condition
Condition<-name


# Create a dataframe with file paths and file names
df<-data.frame(Sample_Name=name,
           path_fwdReads=paste0(pathSTD,"/data/raw_data/",fwd),
           path_revReads=paste0(pathSTD,"/data/raw_data/",rev),
           Condition)

# Write a .csv file with fastq file names and paths
write.csv(df, row.names = F,
          file="data/samples_profiler.csv", quote = F)
