# Name: trimmomatic.sh
# Auth: umar.niazi@kcl.ac.uk.
# Date: 31/08/2016
# Desc: shell script to run trimmomatic on a list of files

#!/bin/bash
#$ -S /bin/sh
#$ -pe smp 1 
#$ -cwd
#$ -N trim-test
#$ -j y
#$ -l mem_free=8G
#$ -l h_rt=02:00:00
# check if empty commandline
if [ "$#" -ne 2 ]
then
    echo "usage:"
    exit 0
fi

module load general/JRE/1.8.0_65
module load bioinformatics/trimmomatic/0.36

#INPUT_samples="$@"

java -jar /opt/apps/bioinformatics/trimmomatic/0.36/trimmomatic-0.36.jar PE -phred33 $1 $2 $1_trim.fq.gz $1_unpaired.fq.gz $2_trim.fq.gz $2_unpaired.fq.gz ILLUMINACLIP:/users/k1625253/brc_scratch/Data/MetaData/trimmomatic_adapters.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


