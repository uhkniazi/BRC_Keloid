# download_data.R
# Auth: Umar Niazi 
# Date: 12/08/2016
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
#### dataset Data.id = 1, GSE7890
###############################################################################
## data loading
# load the data, clean and create factors
setwd('GSE7890')
dir.create('Data_external', showWarnings = F)
setwd('Data_external/')
gse =  getGEOSuppFiles('GSE7890')
untar('GSE7890/GSE7890_RAW.tar')
oData = ReadAffy()
setwd('..')
# normalize the data
x.affy = rma(oData)

# get the geo soft file for the sample annotations
gse =  getGEO('GSE7890', GSEMatrix = T, destdir = 'Data_external/')
oExp = gse$GSE7890_series_matrix.txt.gz

# get the samples from the geo object
dfSamples = pData(oExp)

# col names of the affy expression matrix
cn = colnames(exprs(x.affy))
# remove the last .CEL.gz from the names
cn = gsub('(.+)\\.CEL\\.gz', replacement = '\\1', cn, perl = T)

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

## create grouping factors
fDisorder = as.character(pData(x.affy)[,11])
table(gsub('Disorder: (.+)', '\\1', fDisorder, perl=T))
fDisorder = gsub('Disorder: (.+)', '\\1', fDisorder, perl=T)
fDisorder = factor(fDisorder)
table(fDisorder); levels(fDisorder)
x.affy$fDisorder = fDisorder

fCortisone = as.character(pData(x.affy)[,12])
table(gsub('Condition: (.+)', '\\1', fCortisone, perl=T))
fCortisone = gsub('Condition: (.+)', '\\1', fCortisone, perl=T)
fCortisone = factor(fCortisone)
table(fCortisone); levels(fCortisone)
x.affy$fCortisone = fCortisone

#### initial quality checks using PCA
m = exprs(x.affy)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
## choose appropriate factor
fSamples = x.affy$fDisorder
fSamples = x.affy$fCortisone
fSamples = paste

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
dir.create('Objects')
save(x.affy, file='Objects/normalized.object.rds')

############################ apply combat for batch removal


















######################################

# print samples
as.data.frame(table(oExp$source_name_ch1))

# get the whole blood data
i = grep('treatment', x = oExp$source_name_ch1, ignore.case = F, perl = T)
oExp = oExp[,i]

## data normalization
# normalize and log2 transform the data using lumi
oExp.lumi = lumiT(oExp, 'log2')
# remove any NA data
exprs(oExp.lumi) = na.omit(exprs(oExp.lumi))
fSamples = rep(NA, 21)
f = as.character(oExp$source_name_ch1)
# create factors
i = grep('before treatment', f)
fSamples[i] = '0'

i = grep('2 months', f)
fSamples[i] = '2'

i = grep('12 months', f)
fSamples[i] = '12'

fSamples = factor(fSamples, levels = c('12', '2', '0'))
levels(fSamples)
oExp.lumi$fSamples = fSamples

## data normalization
oExp = lumiN(oExp.lumi, method='rsn')
rm(oExp.lumi)
gc()

# add lumi nuIDs 
oExp = addNuID2lumi(oExp, lib.mapping = 'lumiHumanIDMapping' )
# remove any NA data
exprs(oExp) = na.omit(exprs(oExp))
fSamples = oExp$fSamples
table(fSamples)
levels(fSamples)

### perform DE analysis
mDat = exprs(oExp)
# map the nuID to human symbols
cvSym = getSYMBOL(rownames(mDat), 'lumiHumanAll.db')
# number of genes not annotated
table(is.na(cvSym))
# remove unannotated genes
mDat = mDat[!is.na(cvSym),]

design = model.matrix(~fSamples)
colnames(design) = levels(fSamples)
head(design)

fit = lmFit(mDat, design)
fit = eBayes(fit)

# get annotation
df = select(lumiHumanAll.db, keys = rownames(mDat), columns = c('ENTREZID', 'SYMBOL', 'GENENAME'), keytype = 'PROBEID')
# sanity check
nrow(mDat) == nrow(df)
# add annotation to limma object
fit$genes = df
topTable(fit, adjust='BH')

# look at top tables for each comparison
for (i in 2:length(levels(fSamples))){
  print(paste(levels(fSamples)[i], i))
  print(topTable(fit, coef=i, adjust='BH'))
}

# get the list of genes for each comparison i.e. each coefficient compared to base line
lSigGenes.adj = vector('list', length = length(levels(fSamples))-1)
names(lSigGenes.adj) = levels(fSamples)[2:length(levels(fSamples))]

for (i in 2:length(levels(fSamples))){
  p.adj = p.adjust(fit$p.value[,i], method = 'BH')
  lSigGenes.adj[[i-1]] = names(p.adj)[p.adj < 0.01]
}

#cvSigGenes.adj = unique(cvSigGenes.adj)
sapply(lSigGenes.adj, length)

######### Volcano plots
# plot volcano plots
par(p.old)
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1

for (i in seq_along(n)) {
  dfGenes = topTable(fit, coef = n[i], number = Inf)
  f_plotVolcano(dfGenes, paste(names(n[i])), fc.lim = c(-3, 3))
}

# get the common genes
dfRes = topTable(fit, adjust='BH', number=Inf, p.value=0.1)
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1
cvCommonGenes = NULL
for (i in seq_along(n)) {
  cvCommonGenes = append(cvCommonGenes, lSigGenes.adj[[names(n[i])]])
}
cvCommonGenes = unique(cvCommonGenes)

## export the result for cluster analysis
dfData = t(mDat[cvCommonGenes,])
cn = select(lumiHumanAll.db, keys = colnames(dfData), columns = c('ENTREZID'), keytype = 'PROBEID')
colnames(dfData) = cn$ENTREZID
# remove duplicate probes
f = !duplicated(colnames(dfData))
dfData = dfData[,f]
dfData = data.frame(dfData)
dfData$fSamples = fSamples


dir.create('Test_data', showWarnings = F)

write.csv(dfData, file='Test_data/test_data_GSE19491.csv')