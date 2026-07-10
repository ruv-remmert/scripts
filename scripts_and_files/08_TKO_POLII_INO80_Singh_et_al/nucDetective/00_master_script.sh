#!/usr/bin/env bash

export MY_PWD=$(pwd)

#./01_fetch_annotation.sh $MY_PWD
./02_fetch_fastq.sh $MY_PWD
Rscript ./03_createCSV_Profiler.R
Rscript ./04_create_blacklist.R
./05_runNextflow_Profiler.sh $MY_PWD
Rscript ./06_createCSV_Inspector.R
./07_runNextflow_Inspector.sh $MY_PWD
