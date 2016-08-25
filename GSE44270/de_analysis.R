# Name: de_analysis.R
# Auth: Umar Niazi 
# Date: 12/08/2016
# Desc: analysis of the downloaded data

######################### libraries 
library(annotate)
library(limma)
library(NMF)
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CGraphClust/master/CGraphClust.R'
download(url, 'CGraphClust.R')

# load the required packages
source('CGraphClust.R')
# delete the file after source
unlink('CGraphClust.R')

######################## global variables
p.old = par()

########################

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


######################### end functions

## data loading step
# connect to mysql database to find file locations
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
# another way to get the query, preferred
dfSample = dbGetQuery(db, "select * from MetaFile where idData=4;")
# close connection after getting data
dbDisconnect(db)

## pooled lifecycle stages with fisher p.values
loc = paste0(dfSample[dfSample$id == 5, 'location'], dfSample[dfSample$id == 5, 'name'])

load(loc)
annotation(x.affy)
## load the annotation package
library(org.Hs.eg.db)
setwd('GSE44270')
library(GEOquery)
gse = getGEO('GPL6244', destdir = 'Data_external')
# get annotation from this object
dfAnnotation = gse@dataTable@table
cvGenes = dfAnnotation$gene_assignment
i = grep("(NM)|(NR)_", cvGenes, perl=T)
cvGenes = cvGenes[i]
cvGenes = gsub('.*(NM_\\d+).*', '\\1', cvGenes, perl=T)
cvGenes = gsub('.*(NR_\\d+).*', '\\1', cvGenes, perl=T)
dfAnnotation = dfAnnotation[i,]
dfAnnotation$cvGenes = cvGenes

###################################### grouping
dim(x.affy)
# create a subset of data based on cell type
i = which(x.affy$fCell.type == 'fibroblast')
x.affy = x.affy[,i]
dim(x.affy)
mDat = exprs(x.affy)
# replace probe names by enterez ids
rownames(dfAnnotation) = dfAnnotation$ID
table(rownames(mDat) %in% rownames(dfAnnotation))
f = rownames(mDat) %in% rownames(dfAnnotation)
mDat = mDat[f,]
table(rownames(mDat) %in% rownames(dfAnnotation))
rn = dfAnnotation[rownames(mDat), 'cvGenes']
head(rn)
rownames(mDat) = rn
## find matching enterez ids
dfEnterez = AnnotationDbi::select(org.Hs.eg.db, rn, c('ENTREZID', 'SYMBOL'), keytype = 'REFSEQ')
dim(dfEnterez)
dfEnterez = na.omit(dfEnterez)
dim(dfEnterez)
# get only unique ids
f = !duplicated(dfEnterez$REFSEQ)
dfEnterez = dfEnterez[f,]
dim(dfEnterez)
f = !duplicated(dfEnterez$ENTREZID)
dfEnterez = dfEnterez[f,]
dim(dfEnterez)
rownames(dfEnterez) = dfEnterez$REFSEQ
mDat = mDat[dfEnterez$REFSEQ,]
dim(mDat)
rownames(mDat) = dfEnterez$ENTREZID

# sample annotation
dfSamples = pData(x.affy)

# sanity check
identical(rownames(dfSamples), colnames(mDat))

## create grouping factors
fSamples = dfSamples$fDisorder
table(fSamples)
levels(fSamples)
####### end grouping

################################### main analysis section

design = model.matrix(~fSamples)
colnames(design) = levels(fSamples)
head(design)

# ## correlations - treating people as random effects
# corfit = duplicateCorrelation(x.affy, design, block = x.affy$subject_id)
# 
# fit = lmFit(mDat, design, block=x.affy$subject_id, correlation = corfit$consensus.correlation)

fit = lmFit(mDat, design)
fit = eBayes(fit)

# get enterez gene id
# get annotation
df = select(org.Hs.eg.db, keys = rownames(mDat), columns = c('SYMBOL', 'GENENAME'), keytype ='ENTREZID')
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
  lSigGenes.adj[[i-1]] = names(p.adj)[p.adj < 0.1]
}

#cvSigGenes.adj = unique(cvSigGenes.adj)
sapply(lSigGenes.adj, length)

######### Volcano plots
# plot volcano plots
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1

for (i in seq_along(n)) {
  dfGenes = topTable(fit, coef = n[i], number = Inf)
  f_plotVolcano(dfGenes, paste(names(n[i])), fc.lim = c(-4, 4), p.adj.cut = 0.1)
}

### grouping
dfRes = topTable(fit, adjust='BH', number=Inf, p.value=0.1)
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1
cvCommonGenes = NULL
for (i in seq_along(n)) {
  cvCommonGenes = append(cvCommonGenes, lSigGenes.adj[[names(n[i])]])
}
cvCommonGenes = unique(cvCommonGenes)
mCommonGenes = sapply(seq_along(n), function(x) cvCommonGenes %in% lSigGenes.adj[[names(n[x])]])
rownames(mCommonGenes) = cvCommonGenes
colnames(mCommonGenes) = names(n)

#### analysis by grouping genes
# create groups in the data based on 2^n-1 combinations
mCommonGenes.grp = mCommonGenes
set.seed(123)
dm = dist(mCommonGenes.grp, method='binary')
hc = hclust(dm)

# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mCommonGenes.grp = cbind(mCommonGenes.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mCommonGenes.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)
## choose groups where there are at least  genes
#temp3 = temp2[which(temp2[,ncol(temp2)] >= 50),]

# select groups of choice
## early response genes
# rn = rownames(mCommonGenes.grp[cp %in% c(1),])
# #m1 = mDat[rn,]
# df.rn = select(hgu133plus2.db, keys = rn, columns = c('ENTREZID', 'SYMBOL', 'GENENAME'), keytype = 'PROBEID')
# # write csv to look at gene list
# dir.create('Temp')
# write.csv(df.rn[,-1], file=paste('Temp/', 'normal.vs.keloid.only', '.csv', sep=''))
# 
# rn = rownames(mCommonGenes.grp[cp %in% c(2),])
# #m1 = mDat[rn,]
# df.rn = select(hgu133plus2.db, keys = rn, columns = c('ENTREZID', 'SYMBOL', 'GENENAME'), keytype = 'PROBEID')
# # write csv to look at gene list
# write.csv(df.rn[,-1], file=paste('Temp/', 'normal.normalscar.vs.keloid', '.csv', sep=''))
# 
# rn = rownames(mCommonGenes.grp[cp %in% c(3),])
# #m1 = mDat[rn,]
# df.rn = select(hgu133plus2.db, keys = rn, columns = c('ENTREZID', 'SYMBOL', 'GENENAME'), keytype = 'PROBEID')
# # write csv to look at gene list
# write.csv(df.rn[,-1], file=paste('Temp/', 'normalscar.vs.keloid.only', '.csv', sep=''))

################################################################
## all overexpressed genes if interested in
i = 1:nrow(mCommonGenes)

m1 = as.matrix(mCommonGenes[i,])
m1 = mDat[rownames(m1),]

# if using all data
fGroups = fSamples
colnames(m1) = fGroups
m1 = m1[,order(fGroups)]
fGroups = fGroups[order(fGroups)]

# ignore step if stabalization not required
m1 = t(apply(m1, 1, function(x) f_ivStabilizeData(x, fGroups)))
colnames(m1) = fGroups

# scale across rows
rownames(m1) = dfRes[rownames(m1), 'SYMBOL']
mCounts = t(m1)
mCounts = scale(mCounts)
mCounts = t(mCounts)
# threshhold the values
mCounts[mCounts < -3] = -3
mCounts[mCounts > 3] = 3

# draw the heatmap  color='-RdBu:50'
aheatmap(mCounts, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
         annColors=NA, Colv=NA)


# write results csv files
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1

for (i in seq_along(n)) {
  dfGenes = topTable(fit, coef = n[i], number = Inf)
  dfGenes.2 = dfGenes[lSigGenes.adj[[names(n[i])]],]
  rownames(dfGenes.2) = NULL
  f = paste('Temp/', 'Significant_genes_at_1pcFDR_Combined', names(n[i]), '.xls', sep='')
  dfGenes.2 = dfGenes.2[,c(1, 2, 3, 4, 5, 6, 8, 9)]
  write.csv(dfGenes.2, file=f)
}
# 
# ########### pathway analysis using CGraph library
# # uniprot annotation for data
# cvGenes = rownames(mCommonGenes)
# dfGenes = select(org.Hs.eg.db, keys = cvGenes, columns = c('ENTREZID', 'SYMBOL', 'GENENAME', 'UNIPROT'), keytype = 'ENTREZID')
# dfGenes = na.omit(dfGenes)
# 
# # get reactome data
# url = 'http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt'
# dir.create('Data_external', showWarnings = F)
# csReactomeFile = 'Data_external/UniProt2Reactome_All_Levels.txt'
# # download the reactome file if it doesnt exist
# if (!file.exists(csReactomeFile)) download(url, csReactomeFile)
# # map reactome pathways
# dfReactome = read.csv(csReactomeFile, header = F, stringsAsFactors=F, sep='\t')
# x = gsub('\\w+-\\w+-(\\d+)', replacement = '\\1', x = dfReactome$V2, perl = T)
# dfReactome$V2 = x
# dfReactome.sub = dfReactome[dfReactome$V1 %in% dfGenes$UNIPROT,]
# # get the matching positions for uniprot ids in the reactome table
# i = match(dfReactome.sub$V1, dfGenes$UNIPROT)
# dfReactome.sub$ENTREZID = dfGenes$ENTREZID[i]
# dfGraph = dfReactome.sub[,c('ENTREZID', 'V2')]
# dfGraph = na.omit(dfGraph)
# 
# # get expression data
# mCounts = mDat[unique(dfGenes$ENTREZID),]
# # map the probe names to enterez ids
# i = match(rownames(mCounts), dfGenes$ENTREZID)
# rownames(mCounts) = dfGenes$ENTREZID[i]
# fGroups = fSamples
# #fGroups = as.character(fSamples)
# # select only the groups with significant genes
# # n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1
# # i = which(fSamples %in% levels(fSamples)[c(1, n)])
# # fGroups = factor(fGroups[i], levels = levels(fSamples)[c(1, n)])
# # 
# # # subset the count matrix
# # n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1
# # i = which(fSamples %in% levels(fSamples)[c(1, n)])
# # mCounts = mCounts[,i]
# colnames(mCounts) = fGroups
# mCounts = mCounts[,order(fGroups)]
# fGroups = fGroups[order(fGroups)]
# mCounts = t(mCounts)
# 
# # select genes that have a reactome term attached
# n = unique(dfGraph$ENTREZID)
# mCounts = mCounts[,n]
# print(paste('Total number of genes with Reactome terms', length(n)))
# levels(fGroups)
# 
# # create a correlation matrix to decide cor cutoff
# mCor = cor(mCounts)
# 
# # check distribution 
# hist(mCor, prob=T, main='Correlation of genes', xlab='', family='Arial', breaks=20, xaxt='n')
# axis(1, at = seq(-1, 1, by=0.1), las=2)
# 
# # stabalize the data and check correlation again
# mCounts.bk = mCounts
# # stabalize the data
# mCounts.st = apply(mCounts, 2, function(x) f_ivStabilizeData(x, fGroups))
# rownames(mCounts.st) = fGroups
# 
# # create a correlation matrix
# mCor = cor(mCounts.st)
# # check distribution 
# hist(mCor, prob=T, main='Correlation of genes', xlab='', family='Arial', breaks=20, xaxt='n')
# axis(1, at = seq(-1, 1, by=0.1), las=2)
# 
# # create the graph cluster object
# # using absolute correlation vs actual values lead to different clusters
# oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.6, bSuppressPlots = F)
# 
# ## general graph structure
# set.seed(1)
# plot.final.graph(oGr)
# ecount(getFinalGraph(oGr))
# vcount(getFinalGraph(oGr))
# ## community structure
# # plot the main communities and centrality graphs
# ig = getFinalGraph(oGr)
# par(mar=c(1,1,1,1)+0.1)
# set.seed(1)
# ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T, iSize = 10)
# plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
#      vertex.frame.color=NA, mark.groups=NULL, edge.color='lightgrey')
# set.seed(1)
# ig = getFinalGraph(oGr)
# ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 10)
# plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
#      vertex.frame.color=NA, edge.color='darkgrey')
# 
# # get community sizes
# dfCluster = getClusterMapping(oGr)
# colnames(dfCluster) = c('gene', 'cluster')
# rownames(dfCluster) = dfCluster$gene
# # how many genes in each cluster
# iSizes = sort(table(dfCluster$cluster))
# # remove communities smaller than 5 members or choose a size of your liking
# i = which(iSizes <= 5)
# if (length(i) > 0) {
#   cVertRem = as.character(dfCluster[dfCluster$cluster %in% names(i),'gene'])
#   iVertKeep = which(!(V(getFinalGraph(oGr))$name %in% cVertRem))
#   oGr = CGraphClust.recalibrate(oGr, iVertKeep)
# }
# 
# ## centrality diagnostics
# ## centrality parameters should not be correlated significantly and the location of the central
# ## genes can be visualized
# # look at the graph centrality properties
# set.seed(1)
# ig = plot.centrality.graph(oGr)
# # plot the genes or vertex sizes by fold change
# ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 10)
# set.seed(1)
# plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='darkgrey')
# par(p.old)
# 
# ## the diagnostic plots show the distribution of the centrality parameters
# # these diagnostics plots should be looked at in combination with the centrality graphs
# plot.centrality.diagnostics(oGr)
# 
# # get the centrality parameters
# mCent = mPrintCentralitySummary(oGr)
# 
# ## top vertices based on centrality scores
# ## get a table of top vertices 
# dfTopGenes.cent = dfGetTopVertices(oGr, iQuantile = 0.90)
# rownames(dfTopGenes.cent) = dfTopGenes.cent$VertexID
# # assign metadata annotation to these genes and clusters
# dfCluster = getClusterMapping(oGr)
# colnames(dfCluster) = c('gene', 'cluster')
# rownames(dfCluster) = dfCluster$gene
# df = f_dfGetGeneAnnotation(as.character(dfTopGenes.cent$VertexID))
# dfTopGenes.cent = cbind(dfTopGenes.cent[as.character(df$ENTREZID),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
# dfCluster = dfCluster[as.character(dfTopGenes.cent$VertexID),]
# dfTopGenes.cent = cbind(dfTopGenes.cent, Cluster=dfCluster$cluster)
# 
# ####### NOTE: This section of code is very slow, use only if you need data from genbank
# # loop and get the data from genbank
# n = rep(NA, length=nrow(dfTopGenes.cent))
# names(n) = as.character(dfTopGenes.cent$VertexID)
# for (i in seq_along(n)){
#   n[i] = f_csGetGeneSummaryFromGenbank(names(n)[i])
#   # this wait time is required to stop sending queries to ncbi server very quickly
#   Sys.sleep(time = 3)
# }
# cvSum.2 = as.character(dfTopGenes.cent$VertexID)
# dfTopGenes.cent$Summary = n[cvSum.2]
# ####### Section ends
# 
# dir.create('Results', showWarnings = F)
# write.csv(dfTopGenes.cent, file='Temp/Top_Centrality_Genes_duke.csv')
# 
# ## if we want to look at the expression profiles of the top genes
# # plot a heatmap of these top genes
# library(NMF)
# m1 = mCounts[,as.character(dfTopGenes.cent$VertexID)]
# m1 = scale(m1)
# m1 = t(m1)
# # threshhold the values
# m1[m1 < -3] = -3
# m1[m1 > 3] = 3
# rownames(m1) = as.character(dfTopGenes.cent$SYMBOL)
# # draw the heatmap  color='-RdBu:50'
# aheatmap(m1, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
#          annColors=NA, Colv=NA)
# 
# m1 = mCounts[,as.character(dfTopGenes.cent$VertexID)]
# m1 = apply(m1, 2, f_ivStabilizeData, fGroups)
# rownames(m1) = fGroups
# m1 = scale(m1)
# m1 = t(m1)
# # threshhold the values
# m1[m1 < -3] = -3
# m1[m1 > 3] = 3
# rownames(m1) = as.character(dfTopGenes.cent$SYMBOL)
# # draw the heatmap  color='-RdBu:50'
# aheatmap(m1, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
#          annColors=NA, Colv=NA)
# 
# ## in addition to heatmaps the graphs can be plotted
# # plot a graph of these top genes
# # plot for each contrast i.e. base line vs other level
# lev = levels(fGroups)[-1]
# m = mCounts
# m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
# rownames(m) = rownames(mCounts)
# par(mar=c(1,1,1,1)+0.1)
# for(i in 1:length(lev)){
#   ig = induced_subgraph(getFinalGraph(oGr), vids = as.character(dfTopGenes.cent$VertexID))
#   fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
#   ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=10)
#   n = V(ig)$name
#   lab = f_dfGetGeneAnnotation(n)
#   V(ig)$label = as.character(lab$SYMBOL)
#   set.seed(1)
#   plot(ig, vertex.label.cex=0.2, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', 
#        main=paste(lev[i], 'vs', levels(fGroups)[1]))
#   legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
# }
# 
# 
# ### Looking at the largest clique can be informative in the graph
# # plot the graph with location of the clique highlighted
# set.seed(1)
# ig = plot.graph.clique(oGr)
# ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F)
# par(mar=c(1,1,1,1)+0.1)
# set.seed(1)
# plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey')
# 
# # plot the largest clique at each grouping contrast
# lev = levels(fGroups)[-1]
# m = mCounts
# #m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
# #rownames(m) = rownames(mCounts)
# par(mar=c(1,1,1,1)+0.1)
# for(i in 1:length(lev)){
#   ig = induced_subgraph(getFinalGraph(oGr), vids = unlist(getLargestCliques(oGr)))
#   fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
#   ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=10)
#   n = V(ig)$name
#   lab = f_dfGetGeneAnnotation(n)
#   V(ig)$label = as.character(lab$SYMBOL)
#   set.seed(1)
#   plot(ig, layout=layout_with_fr, main=paste(lev[i], 'vs', levels(fGroups)[1]))
#   legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
# }
# 
# 
# 
# ## instead of looking at individual genes we can look at clusters
# ## we can look at the problem from the other direction and look at clusters instead of genes
# # some sample plots
# # mean expression of groups in every cluster
# par(p.old)
# plot.mean.expressions(oGr, t(mCounts), fGroups, legend.pos = 'bottomleft', main='Total Change in Each Cluster', cex.axis=0.7)
# # only significant clusters
# par(mar=c(7, 3, 2, 2)+0.1)
# plot.significant.expressions(oGr, t(mCounts), fGroups, main='Significant Clusters', lwd=1, bStabalize = T, cex.axis=0.7)
# # principal component plots
# pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)
# par(mar=c(4,2,4,2))
# biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)
# # plot summary heatmaps
# # marginal expression level in each cluster
# plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = F)
# plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = T)
# # plot variance of cluster
# m = getSignificantClusters(oGr, t(mCounts), fGroups)$clusters
# 
# csClust = rownames(m)
# length(csClust)
# pdf('Temp/cluster_variance.pdf')
# par(mfrow=c(1,2))
# boxplot.cluster.variance(oGr, m, fGroups, log=T, iDrawCount = length(csClust), las=2)
# for (i in seq_along(csClust)){
#   temp = t(as.matrix(m[csClust[i],]))
#   rownames(temp) = csClust[i]
#   plot.cluster.variance(oGr, temp, fGroups, log=FALSE);
# }
# dev.off(dev.cur())
# #boxplot.cluster.variance(oGr, m, fGroups, log=T, iDrawCount = length(csClust))
# 
# # print the labels of the clusters from reactome table
# i = which(dfReactome.sub$V2 %in% csClust)
# dfCluster.name = dfReactome.sub[i,c('V2', 'V4')]
# dfCluster.name = dfCluster.name[!duplicated(dfCluster.name$V2),]
# rownames(dfCluster.name) = NULL
# dfCluster.name
# 
# #### plot a graph of clusters 
# dfCluster = getClusterMapping(oGr)
# colnames(dfCluster) = c('gene', 'cluster')
# rownames(dfCluster) = dfCluster$gene
# # how many genes in each cluster
# data.frame(sort(table(dfCluster$cluster)))
# #csClust = rownames(m$clusters)
# csClust = as.character(unique(dfCluster$cluster))
# 
# 
# # plot the graphs at each contrast
# lev = levels(fGroups)[-1]
# m = mCounts
# #m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
# #rownames(m) = rownames(mCounts)
# par(mar=c(1,1,1,1)+0.1)
# for(i in 1:length(lev)){
#   ig = getClusterSubgraph(oGr, csClust)
#   #   # plot the largest compoenent only
#   #   com = components(ig)
#   #   com.lar = which.max(com$csize)
#   #   ig = induced_subgraph(ig, vids = V(ig)[which(com$membership == com.lar)])
#   fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
#   ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=10)
#   n = V(ig)$name
#   lab = f_dfGetGeneAnnotation(n)
#   V(ig)$label = as.character(lab$SYMBOL)
#   set.seed(1)
#   plot(ig, vertex.label.cex=0.14, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey',
#        main=paste(lev[i], 'vs', levels(fGroups)[1]))
#   legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
# }
# 
# df = f_dfGetGeneAnnotation(as.character(dfCluster$gene))
# dfCluster = cbind(dfCluster[as.character(df$ENTREZID),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
# write.csv(dfCluster, file='Temp/ClustersToGeneMapping.xls')
# 
# # save the graph and data objects
# keloid_data = list(graph=oGr, matrix=mCounts, groups=fGroups)
# save(keloid_data, file='Objects/keloid_data.rds')
# 
# # saving graph object to visualize in cytoscape or other graph viewers
# lev = levels(fGroups)[-1]
# m = mCounts
# for(i in 1:length(lev)){
#   ig = getClusterSubgraph(oGr, csClust)
#   fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
#   ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=30)
#   n = V(ig)$name
#   lab = f_dfGetGeneAnnotation(n)
#   V(ig)$label = as.character(lab$SYMBOL)
#   nm = paste('Results/keloid_data', lev[i], 'vs', levels(fGroups)[1], '.graphml', sep='')
#   write.graph(ig, file = nm, format = 'graphml')
# }
# 
# 
# #######################################################################################
# ### selection of plots for various clusters
# #######################################################################################
# 
# 
# # Various plots for one cluster of choice
# csClust = '556833'
# 
# lev = levels(fGroups)[-1]
# m = mCounts
# #m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
# #rownames(m) = rownames(mCounts)
# par(mar=c(1,1,1,1)+0.1)
# for(i in 1:length(lev)){
#   ig = getClusterSubgraph(oGr, csClust)
#   fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
#   ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=60)
#   n = V(ig)$name
#   lab = f_dfGetGeneAnnotation(n)
#   V(ig)$label = as.character(lab$SYMBOL)
#   set.seed(1)
#   plot(ig, vertex.label.cex=0.7, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey',
#        main=paste(lev[i], 'vs', levels(fGroups)[1]))
#   legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
# }
# 
# # heatmap of the genes
# ig.sub = getClusterSubgraph(oGr, csClustLabel = csClust)
# n = f_dfGetGeneAnnotation(V(ig.sub)$name)
# mC = t(mCounts)
# mC = mC[n$ENTREZID,]
# rownames(mC) = n$SYMBOL
# mC = t(scale(t(mC)))
# # threshhold the values
# mC[mC < -3] = -3
# mC[mC > +3] = +3
# # draw the heatmap
# hc = hclust(dist(mC))
# aheatmap(mC, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = hc, annRow=NA, 
#          annColors=NA, Colv=NA)
# 
# 
# # if we want to plot variance of one gene at a time
# n = f_dfGetGeneAnnotation(V(ig.sub)$name)
# mC = t(mCounts)
# mC = mC[n$ENTREZID,]
# rownames(mC) = n$SYMBOL
# rn = rownames(mC)
# length(rn)
# 
# for (i in seq_along(rn)){
#   temp = t(as.matrix(mC[rn[i],]))
#   rownames(temp) = rn[i]
#   plot.cluster.variance(oGr, temp, fGroups, log=FALSE);
# }
# 
# ## some enrichment using string db package
library(STRINGdb)
# find database id for human
temp = get_STRING_species(version = '10')
temp[grep('homo', temp$official_name, ignore.case = T, perl=TRUE),]

dbString = STRINGdb$new(version='10', species=9606, input_directory='')
string.mapped = dbString$map(dfGenes.2, 'ENTREZID', removeUnmappedRows = T)

string.mapped = string.mapped[order(string.mapped$adj.P.Val),]
par(p.old)
hits = string.mapped$STRING_id[1:200]
dbString$plot_network(hits)
pdf('Temp/string.pdf')
dbString$plot_network(hits)
dev.off(dev.cur())

# plot enrichment
dbString$plot_ppi_enrichment(string.mapped$STRING_id, quiet = T)

# category for which to compute the enrichment (i.e. "Process", "Component",
# "Function", "KEGG", "Pfam", "InterPro"). The default category is "Process".
enrichmentProcess = dbString$get_enrichment(string.mapped$STRING_id, category = 'Process', iea=F)

## find interactions for patricular proteins
head(string.mapped)
id = which(string.mapped$ENTREZID == '27242')
string.mapped[id,]
# dbString$get_interactions(string.mapped$STRING_id[id])
ig = dbString$get_graph()
# # this graph can be used to do some analysis
