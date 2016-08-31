# Name: trimmomatic.sh
# Auth: umar.niazi@kcl.ac.uk.
# Date: 31/08/2016
# Desc: shell script to run trimmomatic on a list of files

#!/bin/bash

# check if empty commandline
if [ "$#" -ne 2 ]
then
    echo "usage:"
    exit 0
fi

#INPUT_samples="$@"

echo "java -jar trimmomatic-0.35.jar PE -phred33 $1 $2 $1_trim.fq.gz $1_unpaired.fq.gz $2_trim.fq.gz $2_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

