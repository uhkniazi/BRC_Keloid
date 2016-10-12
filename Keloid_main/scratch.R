# File: scratch.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: space for rough work and testing
# Date: 10/10/2016


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

# create the 2 factors
f1 = gsub('(\\w+)-.+', '\\1', dfSample.title$title)
fTitle = factor(f1)

f2 = gsub('Phenotype (\\w+)', '\\1', dfSample.title$group2)
fCondition = factor(f2, levels=c('Control', 'Keloid'))

fTime = gsub('Timepoint (\\d)', '\\1', dfSample.title$group1)
fTime = factor(fTime)

# factors crossed in relation to condition
table(fTitle, fCondition)
# factors nested in relation to time
table(fTitle, fTime)
# create an interaction factor
fCondition.t = factor(fCondition:fTime)
# the factors i.e samples and condition:time are nested
table(fTitle, fCondition.t)

sf = estimateSizeFactorsForMatrix(mCounts.merge)
mCounts.norm = sweep(mCounts.merge, 2, sf, '/')

## perform DE analysis
x = round(mCounts.norm[1,], 0)

library(car)
library(MASS)

d = fitdistr(x, 'normal')
print(d)
qqPlot(x, 'norm')

d = fitdistr(x, 'gamma')
print(d)
qqPlot(x, 'gamma', shape=d$estimate['shape'], rate=d$estimate['rate'])

getalphabeta.poisson = function(lambda){
  m = mean(lambda)
  v = var(lambda)
  alpha = (m^2)/v
  beta = alpha/m
  return(c(alpha=alpha, beta=beta))
}

d = getalphabeta.poisson(x)
qqPlot(x, 'gamma', shape=d['alpha'], rate=d['beta'])

d = fitdistr(x, 'poisson')
qqPlot(x, 'pois', lambda=d$estimate)

d = fitdistr(x, 'negative binomial')

d = qqPlot(x,'nbinom', size=d$estimate['size'], mu=d$estimate['mu'])

dfData = data.frame(resp = x, patient=fTitle, condition=fCondition, time=fTime, cond.time = fCondition.t)

library(multcomp)
fm01 = glm.nb(resp ~ 0 + cond.time, data=dfData)
#fm02 = glm.nb(resp ~ 1 + cond.time, data=dfData, link=identity)
summary(fm01)
fm01.cont = rbind('cont2 vs cont1' = c(-1, 1, 0, 0), 'kelo1 vs cont1' = c(-1, 0, 1, 0))
summary(glht(fm01, fm01.cont))

library(lme4)

#fm01 = glmer.nb(resp ~ 1 + condition + time + (1 | patient), data=dfData)
fm04 = glmer.nb(resp ~ 0 + cond.time + (1 | patient), data=dfData)

#th = fitdistr((dfData$resp), 'negative binomial')
#fm04 = glmer(resp ~ 1 + cond.time + (1 | patient), data = dfData, family=negative.binomial(2, link=log))
summary(fm04)
pr04 = profile(fm04)
library(lattice)
xyplot(pr04)
confint(pr04)
summary(fm04)
splom(pr04)
qqmath(ranef(fm04, condVar=T))
dotplot(ranef(fm04, condVar=T))

qqnorm(resid(fm04))
qqline(resid(fm04))
qqnorm(ranef(fm04)$patient[,1])
qqline(ranef(fm04)$patient[,1])
plot(fitted(fm04), resid(fm04))
lines(lowess(fitted(fm04), resid(fm04)))
plot(fitted(fm04), x)
lines(lowess(fitted(fm04), x))

dfDesign = data.frame(patient = fTitle, condition=fCondition, time=fTime, cond.time = fCondition.t, row.names = colnames(mCounts.merge))
# 
# # which samples have missing values
# cDropSamples = c('N3', 'N6')
# i = which(dfDesign$patient %in% cDropSamples)
# dfDesign = dfDesign[-i,]
# dfDesign = droplevels.data.frame(dfDesign)
# mCounts.merge = mCounts.merge[,-i]

oDseq = DESeqDataSetFromMatrix(mCounts.merge, dfDesign, design = ~ 1 + condition)
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





