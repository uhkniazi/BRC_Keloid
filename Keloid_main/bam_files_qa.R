# File: bam_files_qa.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: quality checks on the bam files
# Date: 15/09/2016


## set variables and source libraries
gcswd = getwd()
setwd('Keloid_main/')
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CBamQuality/master/CBamQuality.R'
download(url, 'CBamQuality.R')

# load the required packages
source('CBamQuality.R')
# delete the file after source
unlink('CBamQuality.R')

## choose the sequence names over which to calculate coverage
cvSeqnames = paste0('chr', 1:22)

## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'File')
# get the query
q = paste0('select Sample.id as sid, Sample.group1 as timepoints, Sample.title, File.* from Sample, File
where (Sample.idData = 2) AND (Sample.group1 like "%Timepoint%") AND (File.idSample = Sample.id AND File.type like "%duplicates removed%" 
AND File.group1 = "S014")')
dfSample = dbGetQuery(db, q)
# close connection after getting data
dbDisconnect(db)
# remove any whitespace from the names
dfSample$name = gsub(" ", "", dfSample$name, fixed = T)
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)

#### set working directory to appropriate location with bam files
setwd('Data_external/Aligned/S014/')
csFiles = list.files('.', pattern = '*.bam')
# check if these files match the file names in database
table(dfSample$name %in% csFiles)

csFiles = dfSample$name

## perform the analysis one sample at a time
## function to write the qa files
write.qa = function(bfn, sn){
  # open bam file
  bf = BamFile(bfn)
  lBam = lapply(sn, function(x){
    return(CBamScaffold(bf, x))
  })
  cat(paste('done', bfn, '\n'))
  return(lBam)
}

lAllBams = lapply(csFiles, function(x){
  write.qa(x, cvSeqnames)
})
names(lAllBams) = dfSample$title
setwd(gcswd)
n = make.names(paste('CBamScaffold rd S014 rds'))
lAllBams$meta = dfSample
lAllBams$desc = paste('CBamScaffold object from S014 sequencing run with quality 10 duplicates removed', date())
n2 = paste0('~/Data/MetaData/', n)
save(lAllBams, file=n2)

library('RMySQL')
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
df = data.frame(idData=2, name=n, type='rds', location='~/Data/MetaData/', 
                comment='CBamScaffold object from S014 sequencing run with quality 10 duplicates removed')
dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)

### create the plots of interest
getwd()
lAllBams$desc = NULL
lAllBams$meta = NULL
pdf(file='Keloid_main/Results/bam.q10.rd.qa.S014.pdf')
par(mfrow=c(2,2))
## returns the coverage for each scaffold in the bam file in matrix
f1 = function(ob){
  mat = sapply(ob, CBamScaffold.getCoverage)
  n = sapply(ob, CBamScaffold.getSeqname)
  colnames(mat) = n
  return(log(mat+1))
}

## call the function f1
lMat = lapply(lAllBams, f1)
i = seq_along(cvSeqnames)

## reorder the matrix, so that same chromosome/scaffold from all bam files are in one matrix
lMat.ordered = lapply(i, function(x){
  m = sapply(lMat, function(y) y[,x])
})

names(lMat.ordered) = cvSeqnames
iCol = rainbow(length(lAllBams))
## plot the lowess profile 
temp = sapply(cvSeqnames, function(x){
  m = apply(lMat.ordered[[x]], 2, function(y) lowess(y, f=2/10)$y )
  matplot(m, type='l', main=paste('Lowess fit to coverage', x), col=iCol, xlab='Bins', sub='Binned Coverage', ylab='Coverage', 
          lty=1:length(iCol))
})

par(mfrow=c(1,1))
plot.new()
legend('center', legend = names(lAllBams), ncol=3, col = iCol, lty=1:length(iCol), cex=0.6, lwd=2)



## average coverage, read width
f1 = function(ob){
  mat = sapply(ob, CBamScaffold.getReadWidthSample)
  n = sapply(ob, CBamScaffold.getSeqname)
  colnames(mat) = n
  return(mean(colMeans(mat)))
}

## call the function f1
ivMat = sapply(lAllBams, f1)
barplot(ivMat, las=2, main='Sequencing run S014 - Average Read Width', ylab = 'Average Read Width', cex.names =0.6)

# number of reads aligned
f1 = function(ob){
  n = sapply(ob, function(x) length(CBamScaffold.getReadWidth(x)))
  n = sum(n)/1e+6
  return(n)
}

iReadCount = sapply(lAllBams, f1)

barplot(iReadCount, las=2, main='Sequencing run S014 - Read Count', ylab = 'No. of Reads in Millions', cex.names =0.6)

# average binned coverage distribution
f1 = function(ob){
  mat = sapply(ob, getCoverageGammaParam)
  n = sapply(ob, CBamScaffold.getSeqname)
  names(mat) = n
  return(mat)
}

## call the function f1
lMat = lapply(lAllBams, f1)

boxplot(lMat, las=2, main='Sequencing run S014 - Average Binned Coverage', ylab = 'Average', outline=F, xaxt='n')
axis(1, at=1:length(lMat), labels = names(lMat), cex.axis=0.7, las=2)

dev.off(dev.cur())



# bf = BamFile('/run/user/1000/gvfs/sftp:host=10.202.64.28,user=k1625253/users/k1625253/Data/ProjectsData/BRC_Keloid/Aligned/S014/K2-2nd_S014_R1_.fastq.gz_q10_sort_rd.bam')
# si = seqinfo(bf)
# si = as.data.frame(si)
# si = si[order(si$seqlengths, decreasing = T),]
# sn = rownames(si[1:24,])