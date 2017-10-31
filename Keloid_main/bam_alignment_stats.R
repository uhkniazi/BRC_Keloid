# File: bam_alignment_stats.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: number of reads aligned in the bam file using samtools
# Date: 15/12/2016

## set wd 
gcswd = getwd()
setwd('Keloid_main/Data_external/Aligned/')

## for each sequencing run get the list of bam files
setwd('S014/')

lFiles = list.files(pattern = 'bam$')

com = sapply(lFiles, function(x) {
  paste('samtools view -c -F 4', x, sep=' ')
})

al.14 = sapply(com, function(x){
  bamc = system(com, intern = T)
  bamc = signif(as.numeric(bamc)/1e+6, 3)
})

## s032
setwd('../S032/')

lFiles = list.files(pattern = 'bam$')

com = sapply(lFiles, function(x) {
  paste('samtools view -c -F 4', x, sep=' ')
})

al.32 = sapply(com, function(x){
  bamc = system(com, intern = T)
  bamc = signif(as.numeric(bamc)/1e+6, 3)
})

# s021
setwd('../S021/')

lFiles = list.files(pattern = 'bam$')

com = sapply(lFiles, function(x) {
  paste('samtools view -c -F 4', x, sep=' ')
})

al.21 = sapply(com, function(x){
  bamc = system(com, intern = T)
  bamc = signif(as.numeric(bamc)/1e+6, 3)
})



