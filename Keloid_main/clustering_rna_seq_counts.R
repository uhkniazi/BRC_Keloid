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
dfSample = dfSample[dfSample$id == 25,]
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

## get the ercc spike in ids
i = grep('^ERCC', rownames(mCounts))
mCounts.spikein = mCounts[i,]

## stop here as most of matrix is zero 

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

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, normalised', ylim=c(-20, 40), xlim=c(-55, -20))
text(pr.out$x[,1:2], labels = dfSample.names$title, pos = 1, cex=0.6)
legend('topright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, normalised', ylim=c(-40, 25), xlim=c(10, 65))
text(pr.out$x[,1:2], labels = dfSample.names$title, pos = 1, cex=0.6)
legend('topright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)

############## some additional QC checks with distribution of gene expressions
## and sample types
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

# adding some noise to remove zeros
n = dim(mCounts)[1] * dim(mCounts)[2]
mCounts.s = mCounts + rnorm(n)

#### standardize the genes first 
s = apply(mCounts.s, 1, sd)
mCounts.s = sweep(mCounts.s, 1, s, '/')

# get the mean vector and total vector for each sample
ivMean = colMeans(mCounts.s)
ivTotal = colSums(mCounts.s)

subject.id = gsub('(K\\d*|N\\d*)-.*', '\\1', dfSample.names$title, perl=T)
dfData = data.frame(ivMean, ivTotal, batch=dfSample.names$group1, condition = (gsub('Phenotype ', '', dfSample.names$phenotype)), 
                    time=dfSample.names$timepoints, subject.id)
library(lattice)
densityplot(~ ivMean, data=dfData, groups=batch, auto.key=TRUE, main='Average Gene Expression Density in Each Batch',
            xlab='Mean Gene Expression')

densityplot(~ ivMean | batch, data=dfData, groups=condition, auto.key=TRUE, main='Average Gene Expression Density in Each Batch',
            xlab='Mean Gene Expression')

densityplot(~ ivMean | batch, data=dfData, groups=time, auto.key=TRUE, main='Average Gene Expression Density in Each Batch',
            xlab='Mean Gene Expression')

dotplot(ivMean ~ subject.id | time, data=dfData, groups=batch, auto.key=TRUE, main='Average Gene Expression in Each Sample',
       xlab='Mean Gene Expression', pch=20, cex.axis=0.7)

f = paste(dfData$condition, dfData$time)
densityplot(~ ivMean | batch, data=dfData, groups=f, auto.key=TRUE, main='Average Gene Expression Density in Each Batch',
            xlab='Mean Gene Expression')

## for normalized data
# adding some noise to remove zeros
n = dim(mCounts)[1] * dim(mCounts)[2]
# adding noise to avoid negative numbers
mCounts.s = mCounts.norm + rnorm(n)

# get the mean vector and total vector for each sample
ivMean = colMeans(mCounts.s)
ivTotal = colSums(mCounts.s)

subject.id = gsub('(K\\d*|N\\d*)-.*', '\\1', dfSample.names$title, perl=T)
dfData = data.frame(ivMean, ivTotal, batch=dfSample.names$group1, condition = (gsub('Phenotype ', '', dfSample.names$phenotype)), 
                    time=dfSample.names$timepoints, subject.id)
densityplot(~ ivMean , data=dfData, groups=batch, auto.key=TRUE, main='Average Gene Expression Density in Each Batch, Normalised',
            xlab='Mean Gene Expression')

densityplot(~ ivMean | batch, data=dfData, groups=condition, auto.key=TRUE, main='Average Gene Expression Density in Each Batch',
            xlab='Mean Gene Expression')

densityplot(~ ivMean | batch, data=dfData, groups=time, auto.key=TRUE, main='Average Gene Expression Density in Each Batch, Normalised',
            xlab='Mean Gene Expression')

dotplot(ivMean ~ subject.id | time, data=dfData, groups=batch, auto.key=TRUE, main='Average Gene Expression in Each Sample, Normalised',
        xlab='Mean Gene Expression', pch=20, cex.axis=0.7)

f = paste(dfData$condition, dfData$time)
densityplot(~ ivMean | batch, data=dfData, groups=f, auto.key=TRUE, main='Average Gene Expression Density in Each Batch, Normalised',
xlab='Mean Gene Expression')


################# merge replicate samples before normalisation 
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

fMerge = factor(dfSample.names$sid)

mCounts.merge = lapply(levels(fMerge), function(x){
  i = which(fMerge == x)
  m = mCounts[,i]
  rowSums(m)
})

mCounts.merge = do.call(cbind, mCounts.merge)

# get the title for these samples from database
##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# get the query

q = paste0('select Sample.* from Sample
where (Sample.id =', as.numeric(levels(fMerge)), ')')
l = lapply(q, function(x) dbGetQuery(db, x))
dfSample.title = do.call(rbind, l)
dfSample.title$title = gsub(' ', '', dfSample.title$title)
# close connection after getting data
dbDisconnect(db)

colnames(mCounts.merge) = dfSample.title$title

dfDesign = data.frame(condition=factor(dfSample.title$group1), row.names = colnames(mCounts.merge))

oDseq = DESeqDataSetFromMatrix(mCounts.merge, dfDesign, design = ~ condition)
oDseq = DESeq(oDseq)
mCounts.norm = counts(oDseq, normalized=T)

mCounts = mCounts.merge
## for normalized data
# adding some noise to remove zeros
n = dim(mCounts)[1] * dim(mCounts)[2]
# adding noise to avoid negative numbers
mCounts.s = t(mCounts.norm + rnorm(n))

# set scaling to TRUE to scale columns 
pr.out = prcomp(mCounts.s, scale = T)

# set the factor for colours
fSamples = factor(dfSample.title$group1)
fSamples = factor(dfSample.title$group2)
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, normalised')
text(pr.out$x[,1:2], labels = dfSample.title$title, pos = 1, cex=0.6)
legend('topleft', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, normalised', ylim=c(-50, 60), xlim=c(-55, 80))
text(pr.out$x[,1:2], labels = dfSample.title$title, pos = 1, cex=0.6)
legend('topright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, normalised', ylim=c(0, 40), xlim=c(-65, -20))
text(pr.out$x[,1:2], labels = dfSample.title$title, pos = 1, cex=0.6)
legend('topright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, normalised', ylim=c(-55, 25), xlim=c(0, 65))
text(pr.out$x[,1:2], labels = dfSample.title$title, pos = 1, cex=0.6)
legend('topright', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)

############ density plots
# adding some noise to remove zeros
n = dim(mCounts)[1] * dim(mCounts)[2]
mCounts.s = mCounts + rnorm(n)

#### standardize the genes first 
s = apply(mCounts.s, 1, sd)
mCounts.s = sweep(mCounts.s, 1, s, '/')

# get the mean vector and total vector for each sample
ivMean = colMeans(mCounts.s)
ivTotal = colSums(mCounts.s)

subject.id = gsub('(K\\d*|N\\d*)-.*', '\\1', dfSample.title$title, perl=T)
dfData = data.frame(ivMean, ivTotal, condition = (gsub('Phenotype ', '', dfSample.title$group2)), 
                    time=dfSample.title$group1, subject.id)
library(lattice)
densityplot(~ ivMean, data=dfData, groups=condition, auto.key=TRUE, main='Average Gene Expression Density in Each Condition',
            xlab='Mean Gene Expression')

densityplot(~ ivMean, data=dfData, groups=time, auto.key=TRUE, main='Average Gene Expression Density in Each Batch',
            xlab='Mean Gene Expression')

dotplot(ivMean ~ subject.id | time, data=dfData, groups=condition, auto.key=TRUE, main='Average Gene Expression in Each Sample',
        xlab='Mean Gene Expression', pch=20, cex.axis=0.7)

f = paste(dfData$condition, dfData$time)
densityplot(~ ivMean, data=dfData, groups=f, auto.key=TRUE, main='Average Gene Expression Density in Each Batch',
            xlab='Mean Gene Expression')

## for normalized data
# adding some noise to remove zeros
n = dim(mCounts)[1] * dim(mCounts)[2]
# adding noise to avoid negative numbers
mCounts.s = mCounts.norm + rnorm(n)

# get the mean vector and total vector for each sample
ivMean = colMeans(mCounts.s)
ivTotal = colSums(mCounts.s)

subject.id = gsub('(K\\d*|N\\d*)-.*', '\\1', dfSample.title$title, perl=T)
dfData = data.frame(ivMean, ivTotal, condition = (gsub('Phenotype ', '', dfSample.title$group2)), 
                    time=dfSample.title$group1, subject.id)

densityplot(~ ivMean, data=dfData, groups=condition, auto.key=TRUE, main='Average Gene Expression Density in Each Condition, Normalised',
            xlab='Mean Gene Expression')

densityplot(~ ivMean, data=dfData, groups=time, auto.key=TRUE, main='Average Gene Expression Density in Each Batch, Normalised',
            xlab='Mean Gene Expression')

dotplot(ivMean ~ subject.id | time, data=dfData, groups=condition, auto.key=TRUE, main='Average Gene Expression in Each Sample, Normalised',
        xlab='Mean Gene Expression', pch=20, cex.axis=0.7)

xyplot(ivMean ~ time | condition, data=dfData, groups=subject.id, auto.key=FALSE, main='Average Gene Expression in Each Subject, Normalised',
       xlab='Time points', pch=20, cex.axis=0.7, type='o', ylab='Mean Gene Expression')


f = paste(dfData$condition, dfData$time)
densityplot(~ ivMean, data=dfData, groups=f, auto.key=TRUE, main='Average Gene Expression Density in Each Batch, Normalised',
            xlab='Mean Gene Expression', from=18, to=32)





