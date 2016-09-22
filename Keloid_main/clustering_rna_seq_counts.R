# File: clustering_rna_seq_counts.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: cluster the samples based on count table
# Date: 20/09/2016


## set variables and source libraries
gcswd = getwd()
setwd('Keloid_main/')
p.old = par()

######################### functions
# Name: f_Plot3DPCA
# Args: mComp = n X 3 matrix with first 3 components as the 3 column vectors
#       color = colours for the points
#       ... additional arguments to the plot function
# Rets: none
# Desc: takes the first 3 components and plots the data in a 3d plot
f_Plot3DPCA = function(mComp, color, ...) {
  x = mComp[,1]
  y = mComp[,2]
  z = mComp[,3]
  if (!require(scatterplot3d)) stop('scatterplot3d library required')
  scatterplot3d(x, y, z, color, ...)
}

## connect to mysql database to get find path to appropriate file
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# get the query
q = paste0('select MetaFile.* from MetaFile
where (MetaFile.idData = 2) AND (MetaFile.comment like "%count%")')
dfSample = dbGetQuery(db, q)
q = paste0('select Sample.id as sid, Sample.group1 as timepoints, Sample.group2 as phenotype, Sample.title, File.* from Sample, File
where (Sample.idData = 2) AND (Sample.group1 like "%Timepoint%") AND (File.idSample = Sample.id AND File.type like "%duplicates removed%")')
dfSample.names = dbGetQuery(db, q)
# close connection after getting data
dbDisconnect(db)

n = paste0(dfSample$location, dfSample$name)

load(n)

## make count matrix
names(lCounts)
mCounts = do.call(cbind, lCounts)

# reorder the count matrix columns according to the order in samples table
i = match(dfSample.names$name, colnames(mCounts))
mCounts = mCounts[,i]
# sanity check
identical(dfSample.names$name, colnames(mCounts))

#### scale each variable vector i.e. gene
## add a normal jitter to each cell to remove zeros
n = dim(mCounts)[1] * dim(mCounts)[2]
mCounts = mCounts + rnorm(n)

#### standardize samples first
s = apply(mCounts, 2, sd)
mCounts.s = sweep(mCounts, 2, s, '/')

## PCA with strandardizing samples
mCounts.s = t(mCounts.s)
# set scaling to TRUE to scale variables i.e. genes in columns
pr.out = prcomp(mCounts.s, scale = T)

# set the factor for colours
fSamples = factor(dfSample.names$group1)
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, not normalized')
text(pr.out$x[,1:2], labels = dfSample.names$title, pos = 1, cex=0.6)
legend('topright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)

## normalize these datasets with DESeq2 
library(DESeq2)

n = paste0(dfSample$location, dfSample$name)

load(n)

## make count matrix
names(lCounts)
mCounts = do.call(cbind, lCounts)

# reorder the count matrix columns according to the order in samples table
i = match(dfSample.names$name, colnames(mCounts))
mCounts = mCounts[,i]
# sanity check
identical(dfSample.names$name, colnames(mCounts))

dfDesign = data.frame(condition=factor(dfSample.names$timepoints), row.names = colnames(mCounts))

oDseq = DESeqDataSetFromMatrix(mCounts, dfDesign, design = ~ condition)
mCounts.rlog = assays(rlog(oDseq))
oDseq = DESeq(oDseq)
mCounts.norm = counts(oDseq, normalized=T)
mCounts.rlog = mCounts.rlog@listData[[1]]

# adding some noise to remove zeros
n = dim(mCounts)[1] * dim(mCounts)[2]
mCounts.s = t(mCounts.rlog + rnorm(n))

# set scaling to TRUE to scale columns 
pr.out = prcomp(mCounts.s, scale = T)

# set the factor for colours
fSamples = factor(dfSample.names$group1)
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, log transformed')
text(pr.out$x[,1:2], labels = dfSample.names$title, pos = 1, cex=0.6)
legend('topright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)

## for normalized data
# adding some noise to remove zeros
n = dim(mCounts)[1] * dim(mCounts)[2]
# adding noise to avoid negative numbers
mCounts.s = t(mCounts.norm + rnorm(n))

# set scaling to TRUE to scale columns 
pr.out = prcomp(mCounts.s, scale = T)

# set the factor for colours
fSamples = factor(dfSample.names$group1)
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, normalised')
text(pr.out$x[,1:2], labels = dfSample.names$title, pos = 1, cex=0.6)
legend('topright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, normalised', ylim=c(-40, 60), xlim=c(-55, 80))
text(pr.out$x[,1:2], labels = dfSample.names$title, pos = 1, cex=0.6)
legend('topright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)




