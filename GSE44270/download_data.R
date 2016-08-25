# File: download_data.R
# Auth: Umar Niazi 
# Date: 24/08/2016
# Desc: analysis of dataset to generate test data

library(GEOquery)
library(Biobase)
library(affy)
# # library(lumi)
# # library(lumiHumanAll.db)
# # library(lumiHumanIDMapping)
# library(annotate)
# library(limma)

## internal function
# Function: f_lGetPCAClusterCount
# Desc: Takes first two components of the PCA and counts possible clusters. The function does this by 
#       binning the vector of data into X bins, assigning a class label to each bin, counting how many
#       observations in each bin, any total number of bins with at least one observations. This is calculated
#       for both the components of the pca matrix, and the max number of bins with at least one observation, in
#       first or second dimension is reported along with a data.frame with cluster labels.
# Args: pr.out = principal component object returned by prcomp function
# Rets: returns list with 2 elements: 
#       1 - cluster.count = possible number of clusters in the data
#       2 - cluster.label = data.frame with cluster labels
f_lGetPCAClusterCount = function(pr.out){
  # how many clusters in data, using first 2 components
  x1 = pr.out$x[,1]
  x2 = pr.out$x[,2]
  # bin the data from the 2 components
  h1 = hist(x1, plot=F)
  # give a class label to each bin
  c1 = cut(x1, h1$breaks, labels = 1:(length(h1$mids)))
  h2 = hist(x2, plot=F)
  c2 = cut(x2, h2$breaks, labels = 1:(length(h2$mids)))
  # labels for vectors and the class labels
  dfClust = data.frame(lab=names(x1), c1, c2)
  # get contingency table
  mClust = as.matrix(table(c1 = dfClust$c1, c2 = dfClust$c2))
  # count the max of row and col sums that are not zero
  ir = length(which(rowSums(mClust) != 0))
  ic = length(which(colSums(mClust) != 0))
  iClust.count = ifelse(ir > ic, ir, ic)
  lRet = list(cluster.count=iClust.count, cluster.label=dfClust)
  return(lRet)
}

f_plotVolcano = function(dfGenes, main, p.adj.cut = 0.1, fc.lim = c(-3, 3)){
  p.val = -1 * log10(dfGenes$P.Value)
  fc = dfGenes$logFC
  # cutoff for p.value y.axis
  y.cut = -1 * log10(0.01)
  col = rep('lightgrey', times=length(p.val))
  c = which(dfGenes$adj.P.Val < p.adj.cut)
  col[c] = 'red'
  plot(fc, p.val, pch=20, xlab='Fold Change', ylab='-log10 P.Value', col=col, main=main, xlim=fc.lim)
  abline(v = 0, col='grey', lty=2)
  abline(h = y.cut, col='red', lty=2)
  # second cutoff for adjusted p-values
  y.cut = quantile(p.val[c], probs=0.95)
  abline(h = y.cut, col='red')
  # identify these genes
  g = which(p.val > y.cut)
  lab = dfGenes[g, 'SYMBOL']
  text(dfGenes$logFC[g], y = p.val[g], labels = lab, pos=2, cex=0.6)
}
### end internal functions

# global variables
p.old = par()

###############################################################################
#### dataset Data.id = 4, GSE44270
###############################################################################
## data loading
# load the data, clean and create factors
setwd('GSE44270')
dir.create('Data_external', showWarnings = F)
setwd('Data_external/')
gse =  getGEOSuppFiles('GSE44270')
untar('GSE44270/GSE44270_RAW.tar')
oData = ReadAffy()
setwd('..')
# normalize the data
x.affy = rma(oData)

# get the geo soft file for the sample annotations
gse.soft =  getGEO('GSE44270', GSEMatrix = T, destdir = 'Data_external/')
oExp = gse.soft$GSE44270_series_matrix.txt.gz

# get the samples from the geo object
dfSamples = pData(oExp)

# col names of the affy expression matrix
cn = colnames(exprs(x.affy))
# remove the last .CEL.gz from the names
cn = gsub('(GSM\\d+)_.+\\.CEL\\.gz', replacement = '\\1', cn, perl = T)

# order the sample names according to cel files
table(cn %in% rownames(dfSamples))
i = match(cn, rownames(dfSamples))
dfSamples = dfSamples[i,]

## complete the expression set object
pData(x.affy) = dfSamples
colnames(exprs(x.affy)) = cn
# sanity check
t1 = colnames(x.affy)
t2 = rownames(pData(x.affy))
table(t1 %in% t2)
all(t1 == t2)
## create grouping factors
cvSource = as.character(pData(x.affy)[,'source_name_ch1'])
table(cvSource)
## create factors for cell type
fCell.type = rep(NA, length=length(cvSource))
i = grep('fibroblast', cvSource, ignore.case = T)
fCell.type[i] = 'fibroblast'
fCell.type[-i] = 'keratinocyte'
fCell.type = factor(fCell.type)
## keloid vs normal individual
fDisorder = rep(NA, length=length(cvSource))
i = grep('keloid', cvSource, ignore.case = T)
fDisorder[i] = 'keloid'
i = grep('without keloid', cvSource, ignore.case = T)
fDisorder[i] = 'healthy'
i = grep('non-lesion', cvSource, ignore.case = T)
fDisorder[i] = 'keloid.non-lesion'
# sanity check
data.frame(fDisorder, cvSource)
fDisorder = factor(fDisorder, levels = c('keloid', 'keloid.non-lesion', 'healthy'))
table(fDisorder); levels(fDisorder)
x.affy$fDisorder = fDisorder
x.affy$fCell.type = fCell.type

#### initial quality checks using PCA
m = exprs(x.affy)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
## choose appropriate factor
fSamples = x.affy$fDisorder
fSamples = x.affy$fCell.type

col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]
# plot the pca components
par(mfrow=c(2,2))
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], ncol = 2)
plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, normalized')
plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3, normalized')
plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3, normalized')
par(p.old)

getwd()
dir()
n = make.names(paste('GSE44270', 'affymetrix expression set normalized object.rds'))
n2 = paste0('Objects/', n)
save(x.affy, file=n2)

## access the mysql database to create entry for this data object
library('RMySQL')
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
df = data.frame(idData=4, name=n, type='rds', location='~/Data/MetaData/', comment='affymetrix normalised expression data with addded factors')
dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)

