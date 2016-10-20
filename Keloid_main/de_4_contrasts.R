# File: de_4_contrasts.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 20/10/2016
# Desc: DE analysis for the keloid data set

## set variables and source libraries
gcswd = getwd()
setwd('Keloid_main/')
p.old = par()

## libraries to load
library('RMySQL')

### space for internal functions

### end of functions


##### connect to mysql database to get count matrix
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# get the query
q = paste0('select MetaFile.* from MetaFile
where (MetaFile.idData = 2) AND (MetaFile.comment like "%count%")')
dfCountMatrix = dbGetQuery(db, q)
# load the sample meta data
q = paste0('select Sample.id as sid, Sample.group1 as timepoints, Sample.group2 as phenotype, Sample.title, File.* from Sample, File
where (Sample.idData = 2) AND (Sample.group1 like "%Timepoint%") AND (File.idSample = Sample.id AND File.type like "%duplicates removed%")')
dfSample.names = dbGetQuery(db, q)
# close connection after getting data
dbDisconnect(db)

################# merge replicate samples before normalisation 
## normalize these datasets with DESeq2 
library(DESeq2)
n = paste0(dfCountMatrix$location, dfCountMatrix$name)
load(n)

## make count matrix
names(lCounts)
mCounts = do.call(cbind, lCounts)

# reorder the count matrix columns according to the order in samples table
i = match(dfSample.names$name, colnames(mCounts))
mCounts = mCounts[,i]
# sanity check
identical(dfSample.names$name, colnames(mCounts))

# merge the duplicate samples i.e. technical replicates
fMerge = factor(dfSample.names$sid)

mCounts.merge = lapply(levels(fMerge), function(x){
  i = which(fMerge == x)
  m = mCounts[,i]
  rowSums(m)
})

mCounts.merge = do.call(cbind, mCounts.merge)

### get the title for these samples from database
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# get the query
q = paste0('select Sample.* from Sample
where (Sample.id =', as.numeric(levels(fMerge)), ')')
l = lapply(q, function(x) dbGetQuery(db, x))
dfSample.title = do.call(rbind, l)
# remove any whitespace from the title
dfSample.title$title = gsub(' ', '', dfSample.title$title)
# close connection after getting data
dbDisconnect(db)

colnames(mCounts.merge) = dfSample.title$title

### factors for comparisons
# title represents the repeated samples from individuals
f1 = gsub('(\\w+)-.+', '\\1', dfSample.title$title)
fTitle = factor(f1)

# condition represents the biological condition
f2 = gsub('Phenotype (\\w+)', '\\1', dfSample.title$group2)
fCondition = factor(f2, levels=c('Control', 'Keloid'))

# time represents the timing of sample collection
fTime = gsub('Timepoint (\\d)', '\\1', dfSample.title$group1)
fTime = factor(fTime)

## how are the factors related to each other i.e. crossed or nested
# factors are nested in relation to condition title and condition
# each subject has only one condition. this is important as it means deseq will not handle this data
# unless these factors were crossed
table(fTitle, fCondition)
# factors are crossed in relation to condition and time
# each subject was measured at 2 time points, this can be handeled by deseq
table(fTitle, fTime)
# create an interaction factor between time and condition
fCondition.t = factor(fCondition:fTime)
# the factors are partially nested and partially crossed?
table(fTitle, fCondition.t)

dfSample.title$fTitle = fTitle
dfSample.title$fCondition.t = fCondition.t
str(dfSample.title)
# put all the data - matrix and annotation in expression set object
rownames(dfSample.title) = dfSample.title$title
identical(rownames(dfSample.title), colnames(mCounts.merge))

oExp = ExpressionSet(mCounts.merge)
pData(oExp) = dfSample.title
dim(oExp)

# remove the outlier sample identified in previous analysis, N3-2nd
i = which(oExp$title == 'N3-2nd')
oExp = oExp[,-i]
dim(oExp)

# normalise the count matrix using deseq method
sf = estimateSizeFactorsForMatrix(exprs(oExp))
exprs(oExp) = sweep(exprs(oExp), 2, sf, '/')
oExp$LibrarySizeFactor = sf
# drop any unused levels
pData(oExp) = droplevels.data.frame(pData(oExp))

## save this object in the database for future use
## NOTE: don't run this segment of code again as object is already saved
## commenting for safety
# n = make.names(paste('Expression set object for keloids q10 rdup rds'))
# n2 = paste0('~/Data/MetaData/', n)
# save(oExp, file=n2)
# 
# library('RMySQL')
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=2, name=n, type='rds', location='~/Data/MetaData/', 
#                 comment='Normalised Expression set object with comparison factor list for keloids S014 S021 and S032 sequencing runs with quality 10 duplicates removed')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)


