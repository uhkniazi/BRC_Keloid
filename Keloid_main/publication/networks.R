# File: networks.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 16/11/2017
# Desc: generating some figures for network analysis

## set variables and source libraries
gcswd = getwd()
setwd('Keloid_main/publication/')
p.old = par()

## libraries to load
library('RMySQL')

##### connect to mysql database to get count matrix
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# get the query
q = paste0('select MetaFile.* from MetaFile
where (MetaFile.idData = 2)')
dfMetaFiles = dbGetQuery(db, q)
## load the required data
n = paste0(dfMetaFiles$location, dfMetaFiles$name)
i = which(dfMetaFiles$id %in% c(21, 22, 23))
n = n[i]
for (i in seq_along(n)) load(n[i])

# close connection after getting data
dbDisconnect(db)

library(multcomp)
library(glmmADMB)

lGlm.sub = lGlm[!sapply(lGlm, is.null)]

## remove any null elements
table(sapply(lGlm.rep, is.null))
f = sapply(lGlm.rep, is.null)
lGlm.rep[f] = NULL

i = match(names(lGlm.rep), names(lGlm.sub))
head(names(lGlm.sub[i]),3)
head(names(lGlm.rep), 3)
lGlm.sub[i] = lGlm.rep

##################################################################
### repeat the contrasts 
# extract a contrast at a time
mContrasts = rbind('Control:2 vs Control:1' = c(-1, 1, 0, 0),
                   'Keloid:1 vs Control:1' = c(-1, 0, 1, 0),
                   #'Keloid:2 vs Control:1' = c(-1, 0, 0, 1),
                   'Keloid:2 vs Keloid:1' = c(0, 0, -1, 1),
                   'Keloid:2 vs Control:2' = c(0, -1, 0, 1))


## perform contrasts tests
index = 1:length(lGlm.sub)

lContrast1 = lapply(index, function(dat){
  tryCatch({
    s = summary(glht(lGlm.sub[[dat]], t(mContrasts[1,])))
    ret = c(s$test$coefficients[1], s$test$pvalues[1])
    names(ret) = c('logfc', 'p.value')
    return(ret)
  }, error=function(e) {
    ret = c('logfc'=NA, 'p.value'=NA)
    return(ret)
  })
})

dfContrast1 = data.frame(do.call(rbind, lContrast1))
dfContrast1$p.adj = p.adjust(dfContrast1$p.value, method = 'BH')
rownames(dfContrast1) = names(lGlm.sub)

# second contrast
lContrast2 = lapply(index, function(dat){
  tryCatch({
    s = summary(glht(lGlm.sub[[dat]], t(mContrasts[2,])))
    ret = c(s$test$coefficients[1], s$test$pvalues[1])
    names(ret) = c('logfc', 'p.value')
    return(ret)
  }, error=function(e) {
    ret = c('logfc'=NA, 'p.value'=NA)
    return(ret)
  })
})

dfContrast2 = data.frame(do.call(rbind, lContrast2))
dfContrast2$p.adj = p.adjust(dfContrast2$p.value, method = 'BH')
rownames(dfContrast2) = names(lGlm.sub)

# third contrast
lContrast3 = lapply(index, function(dat){
  tryCatch({
    s = summary(glht(lGlm.sub[[dat]], t(mContrasts[3,])))
    ret = c(s$test$coefficients[1], s$test$pvalues[1])
    names(ret) = c('logfc', 'p.value')
    return(ret)
  }, error=function(e) {
    ret = c('logfc'=NA, 'p.value'=NA)
    return(ret)
  })
})

dfContrast3 = data.frame(do.call(rbind, lContrast3))
dfContrast3$p.adj = p.adjust(dfContrast3$p.value, method = 'BH')
rownames(dfContrast3) = names(lGlm.sub)

# fourth contrast
lContrast4 = lapply(index, function(dat){
  tryCatch({
    s = summary(glht(lGlm.sub[[dat]], t(mContrasts[4,])))
    ret = c(s$test$coefficients[1], s$test$pvalues[1])
    names(ret) = c('logfc', 'p.value')
    return(ret)
  }, error=function(e) {
    ret = c('logfc'=NA, 'p.value'=NA)
    return(ret)
  })
})

dfContrast4 = data.frame(do.call(rbind, lContrast4))
dfContrast4$p.adj = p.adjust(dfContrast4$p.value, method = 'BH')
rownames(dfContrast4) = names(lGlm.sub)

## assign annotation to genes
library(org.Hs.eg.db)

df = select(org.Hs.eg.db, as.character(rownames(dfContrast1)), c('SYMBOL', 'GENENAME'), 'ENTREZID')
table(duplicated(df$ENTREZID))
dfContrast1$SYMBOL = df$SYMBOL
dfContrast1$GENENAME = df$GENENAME

df = select(org.Hs.eg.db, as.character(rownames(dfContrast2)), c('SYMBOL', 'GENENAME'), 'ENTREZID')
table(duplicated(df$ENTREZID))
dfContrast2$SYMBOL = df$SYMBOL
dfContrast2$GENENAME = df$GENENAME

df = select(org.Hs.eg.db, as.character(rownames(dfContrast3)), c('SYMBOL', 'GENENAME'), 'ENTREZID')
table(duplicated(df$ENTREZID))
dfContrast3$SYMBOL = df$SYMBOL
dfContrast3$GENENAME = df$GENENAME

df = select(org.Hs.eg.db, as.character(rownames(dfContrast4)), c('SYMBOL', 'GENENAME'), 'ENTREZID')
table(duplicated(df$ENTREZID))
dfContrast4$SYMBOL = df$SYMBOL
dfContrast4$GENENAME = df$GENENAME

## estimate dispersions and some quality plots
par(mfrow=c(2,2))
calculateDispersion = function(fm){
  n = length(resid(fm))
  sqrt(sum(c(as.numeric(resid(fm)), as.numeric(fm$U[[1]]))^2)/n)
}

iDispersion = sapply(lGlm.sub, calculateDispersion)

plotDispersion = function(dis, dat, p.cut, title){
  col = rep('grey', length.out=(nrow(dat)))
  col[dat$p.adj < p.cut] = 'red'
  plot(dis, dat$logfc, col=col, pch=20, main=title, xlab='Dispersion', ylab='logFC')
}

plotDispersion(iDispersion, dfContrast1, 0.01, 'Control:2 vs Control:1')
plotDispersion(iDispersion, dfContrast2, 0.1, 'Keloid:1 vs Control:1')
plotDispersion(iDispersion, dfContrast3, 0.01, 'Keloid:2 vs Keloid:1')
plotDispersion(iDispersion, dfContrast4, 0.1, 'Keloid:2 vs Control:2')

iMean = rowMeans(exprs(oExp)[names(lGlm.sub),])

plotMeanFC = function(m, dat, p.cut, title){
  col = rep('grey', length.out=(nrow(dat)))
  col[dat$p.adj < p.cut] = 'red'
  #mh = cut(m, breaks=quantile(m, 0:50/50), include.lowest = T)
  plot(m, dat$logfc, col=col, pch=20, main=title, xlab='log Mean', ylab='logFC')
}

par(mfrow=c(2,2))
plotMeanFC(log(iMean), dfContrast1, 0.01, 'Control:2 vs Control:1')
plotMeanFC(log(iMean), dfContrast2, 0.1, 'Keloid:1 vs Control:1')
plotMeanFC(log(iMean), dfContrast3, 0.01, 'Keloid:2 vs Keloid:1')
plotMeanFC(log(iMean), dfContrast4, 0.1, 'Keloid:2 vs Control:2')

# add dispersion parameter to ecah gene
dfContrast1$Dispersion = iDispersion[rownames(dfContrast1)]
dfContrast2$Dispersion = iDispersion[rownames(dfContrast2)]
dfContrast3$Dispersion = iDispersion[rownames(dfContrast3)]
dfContrast4$Dispersion = iDispersion[rownames(dfContrast4)]

### grouping of genes
dfContrast1.sub = na.omit(dfContrast1[dfContrast1$p.adj < 0.05 & dfContrast1$Dispersion > 0.4,])
dfContrast2.sub = na.omit(dfContrast2[dfContrast2$p.adj < 0.05 & dfContrast2$Dispersion > 0.4,])
dfContrast3.sub = na.omit(dfContrast3[dfContrast3$p.adj < 0.05 & dfContrast3$Dispersion > 0.4,])
dfContrast4.sub = na.omit(dfContrast4[dfContrast4$p.adj < 0.05 & dfContrast4$Dispersion > 0.4,])

#cvCommonGenes = unique(c(rownames(dfContrast1.sub), rownames(dfContrast2.sub), rownames(dfContrast3.sub), rownames(dfContrast4.sub)))
# use only genes from the contrast of interest
cvCommonGenes = rownames(dfContrast3.sub)

## download the cgraph library
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CGraphClust/master/CGraphClust.R'
download(url, 'CGraphClust.R')

# load the required packages
source('CGraphClust.R')
# delete the file after source
unlink('CGraphClust.R')

####################################################################################################
########### pathway analysis using CGraph library
# uniprot annotation for data
cvGenes = cvCommonGenes
dfGenes = select(org.Hs.eg.db, keys = cvGenes, columns = c('ENTREZID', 'SYMBOL', 'GENENAME', 'UNIPROT'), keytype = 'ENTREZID')
dfGenes = na.omit(dfGenes)

# get reactome data
url = 'http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt'
dir.create('Data_external', showWarnings = F)
csReactomeFile = 'Data_external/UniProt2Reactome_All_Levels.txt'
# download the reactome file if it doesnt exist
if (!file.exists(csReactomeFile)) download(url, csReactomeFile)
# map reactome pathways
dfReactome = read.csv(csReactomeFile, header = F, stringsAsFactors=F, sep='\t')
x = gsub('.+-(\\d+)', replacement = '\\1', x = dfReactome$V2, perl = T)
dfReactome$V2 = x
dfReactome.sub = dfReactome[dfReactome$V1 %in% dfGenes$UNIPROT,]
# get the matching positions for uniprot ids in the reactome table
i = match(dfReactome.sub$V1, dfGenes$UNIPROT)
dfReactome.sub$ENTREZID = dfGenes$ENTREZID[i]
dfGraph = dfReactome.sub[,c('ENTREZID', 'V2')]
dfGraph = na.omit(dfGraph)
head(dfGraph)
# get expression data after subsetting for contrast of interest
i = which((oExp)$fCondition.t %in% c('Keloid:1', 'Keloid:2'))
oExp.sub = oExp[,i]
pData(oExp.sub) = droplevels.data.frame(pData(oExp.sub))
mCounts = exprs(oExp.sub)[unique(dfGenes$ENTREZID),]
dim(mCounts)

fGroups = oExp.sub$fCondition.t
names(fGroups) = oExp.sub$fTitle
colnames(mCounts) = fGroups
# reorder on grouping factor
mCounts = mCounts[,order(fGroups)]
fGroups = fGroups[order(fGroups)]
mCounts = t(mCounts)

# select genes that have a reactome term attached
n = unique(dfGraph$ENTREZID)
mCounts = mCounts[,n]
dim(mCounts)
levels(fGroups)

# create a correlation matrix to decide cor cutoff
mCor = cor(log(mCounts+0.5))

# check distribution 
hist(mCor, prob=T, main='Correlation of genes', xlab='', family='Arial', breaks=20, xaxt='n')
axis(1, at = seq(-1, 1, by=0.1), las=2)

# create the graph cluster object
# using absolute correlation vs actual values lead to different clusters
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.6, bSuppressPlots = F)

## general graph structure
set.seed(1)
plot.final.graph(oGr)
ecount(getFinalGraph(oGr))
vcount(getFinalGraph(oGr))
## community structure
# plot the main communities and centrality graphs
ig = getFinalGraph(oGr)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T, iSize = 2)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, mark.groups=NULL, edge.color='lightgrey')
set.seed(1)
ig = getFinalGraph(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 2)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, edge.color='darkgrey')

# get community sizes
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
dim(dfCluster)
str(dfCluster)
# how many genes in each cluster
iSizes = sort(table(dfCluster$cluster))
# remove communities smaller than 5 members or choose a size of your liking
i = which(iSizes <= 5)
if (length(i) > 0) {
  cVertRem = as.character(dfCluster[dfCluster$cluster %in% names(i),'gene'])
  iVertKeep = which(!(V(getFinalGraph(oGr))$name %in% cVertRem))
  oGr = CGraphClust.recalibrate(oGr, iVertKeep)
}

## centrality diagnostics
## centrality parameters should not be correlated significantly and the location of the central
## genes can be visualized
# look at the graph centrality properties
set.seed(1)
ig = plot.centrality.graph(oGr)
# plot the genes or vertex sizes by fold change
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 2)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='darkgrey')
par(p.old)

## the diagnostic plots show the distribution of the centrality parameters
# these diagnostics plots should be looked at in combination with the centrality graphs
plot.centrality.diagnostics(oGr)

# get the centrality parameters
mCent = mPrintCentralitySummary(oGr)

## top vertices based on centrality scores
## get a table of top vertices 
dfTopGenes.cent = dfGetTopVertices(oGr, iQuantile = 0.90)
rownames(dfTopGenes.cent) = dfTopGenes.cent$VertexID
dim(dfTopGenes.cent)
# assign metadata annotation to these genes and clusters
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
df = f_dfGetGeneAnnotation(as.character(dfTopGenes.cent$VertexID))
dfTopGenes.cent = cbind(dfTopGenes.cent[as.character(df$ENTREZID),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
dfCluster = dfCluster[as.character(dfTopGenes.cent$VertexID),]
dfTopGenes.cent = cbind(dfTopGenes.cent, Cluster=dfCluster$cluster)

####### NOTE: This section of code is very slow, use only if you need data from genbank
# loop and get the data from genbank
n = rep(NA, length=nrow(dfTopGenes.cent))
names(n) = as.character(dfTopGenes.cent$VertexID)
for (i in seq_along(n)){
  n[i] = f_csGetGeneSummaryFromGenbank(names(n)[i])
  # this wait time is required to stop sending queries to ncbi server very quickly
  Sys.sleep(time = 3)
}
cvSum.2 = as.character(dfTopGenes.cent$VertexID)
dfTopGenes.cent$Summary = n[cvSum.2]
####### Section ends

write.csv(dfTopGenes.cent, file='results/Top_Centrality_Genes.xls')

## if we want to look at the expression profiles of the top genes
# plot a heatmap of these top genes
library(NMF)
m1 = log(mCounts[,as.character(dfTopGenes.cent$VertexID)]+0.5)
m1 = scale(m1)
m1 = t(m1)
# threshhold the values
m1[m1 < -3] = -3
m1[m1 > 3] = 3
rownames(m1) = as.character(dfTopGenes.cent$SYMBOL)
# draw the heatmap  color='-RdBu:50'
aheatmap(m1, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
         annColors=NA, Colv=NA)

m1 = log(mCounts[,as.character(dfTopGenes.cent$VertexID)]+0.5)
m1 = apply(m1, 2, f_ivStabilizeData, fGroups)
rownames(m1) = fGroups
m1 = scale(m1)
m1 = t(m1)
# threshhold the values
m1[m1 < -3] = -3
m1[m1 > 3] = 3
rownames(m1) = as.character(dfTopGenes.cent$SYMBOL)
# draw the heatmap  color='-RdBu:50'
aheatmap(m1, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
         annColors=NA, Colv=NA)

## in addition to heatmaps the graphs can be plotted
# plot a graph of these top genes
# plot for each contrast i.e. base line vs other level
lev = levels(fGroups)[-1]
m = mCounts
m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
rownames(m) = rownames(mCounts)
par(mar=c(1,1,1,1)+0.1)
for(i in 1:length(lev)){
  ig = induced_subgraph(getFinalGraph(oGr), vids = as.character(dfTopGenes.cent$VertexID))
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=5)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  set.seed(1)
  plot(ig, vertex.label.cex=0.2, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', 
       main=paste(lev[i], 'vs', levels(fGroups)[1]))
  legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
}


### Looking at the largest clique can be informative in the graph
# plot the graph with location of the clique highlighted
set.seed(1)
ig = plot.graph.clique(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 2)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey')

# plot the largest clique at each grouping contrast
lev = levels(fGroups)[-1]
m = mCounts
#m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
#rownames(m) = rownames(mCounts)
par(mar=c(1,1,1,1)+0.1)
for(i in 1:length(lev)){
  ig = induced_subgraph(getFinalGraph(oGr), vids = unlist(getLargestCliques(oGr)))
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=2)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  set.seed(1)
  plot(ig, layout=layout_with_fr, main=paste(lev[i], 'vs', levels(fGroups)[1]))
  legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
}

## instead of looking at individual genes we can look at clusters
## we can look at the problem from the other direction and look at clusters instead of genes
# some sample plots
# mean expression of groups in every cluster
par(p.old)
plot.mean.expressions(oGr, log(t(mCounts)+0.5), fGroups, legend.pos = 'bottomleft', main='Total Change in Each Cluster', cex.axis=0.7)
# only significant clusters
par(mar=c(7, 3, 2, 2)+0.1)
plot.significant.expressions(oGr, log(t(mCounts)+0.5), fGroups, main='Significant Clusters', lwd=1, bStabalize = T, cex.axis=0.7)#, p.cut=0.05 )
# principal component plots
pr.out = plot.components(oGr, log(t(mCounts)+0.5), fGroups, bStabalize = T)#, p.cut=0.05)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)
# plot summary heatmaps
# marginal expression level in each cluster
plot.heatmap.significant.clusters(oGr, log(t(mCounts)+0.5), fGroups, bStabalize = F)#, p.cut=0.05)
#plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = T, p.cut=0.05)
# plot variance of cluster
#m = getSignificantClusters(oGr, log(t(mCounts)+0.5), fGroups)$clusters
m = getClusterMarginal(oGr, log(t(mCounts)+0.5))
csClust = rownames(m)
length(csClust)
#pdf('Temp/cluster_variance.pdf')
#par(mfrow=c(1,2))
boxplot.cluster.variance(oGr, m, fGroups, log=T, iDrawCount = length(csClust), las=2)
# for (i in seq_along(csClust)){
#   temp = t(as.matrix(m[csClust[i],]))
#   rownames(temp) = csClust[i]
#   plot.cluster.variance(oGr, temp, fGroups, log=FALSE);
# }
# dev.off(dev.cur())


# print the labels of the clusters from reactome table
i = which(dfReactome.sub$V2 %in% csClust)
dfCluster.name = dfReactome.sub[i,c('V2', 'V4')]
dfCluster.name = dfCluster.name[!duplicated(dfCluster.name$V2),]
rownames(dfCluster.name) = NULL
dfCluster.name

#### save a table of clusters to gene mapping
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
write.csv(dfCluster, file='results/dfCluster_members.xls')
# how many genes in each cluster
data.frame(sort(table(dfCluster$cluster)))

mMarginal = getClusterMarginal(oGr, log(t(mCounts)+0.5))

### create lattice plots for each cluster and group i.e. keloids and normals/controls
library(lattice)
## convert to z-scale for each cluster - i.e. vector
dfData = data.frame(scale(t(mMarginal)))
#colnames(dfData) = gsub('X', '', colnames(dfData))
dfData$groups = fGroups
dfData$patient = factor(names(fGroups))
# extract patient conditions from the fGroups
fTime = factor(gsub('^Control:|Keloid:', '', as.character(fGroups)))
fCondition = factor(gsub(':1|:2', '', as.character(fGroups)))
dfData$time = fTime
dfData$condition = fCondition
head(dfData)
str(dfData)
## stack the dataframe expression values for plotting
dfStack = stack(dfData)
str(dfStack)
dfStack$time = dfData$time
dfStack$condition = dfData$condition
dfStack$patient = dfData$patient

xyplot(values ~ time | ind, data=dfStack, type=c('n','smooth'), groups=condition, auto.key=T)
xyplot(values ~ time | ind, data=dfStack, type=c('g', 'p', 'r'), groups=condition, auto.key=T)

## xyplots with population average
mMarginal = getClusterMarginal(oGr, log(t(mCounts)+0.5))
mMarginal = apply(mMarginal, 1, function(x) tapply((x), fGroups, mean))
mMarginal = scale(mMarginal)
dfData = data.frame(mMarginal)
cn = colnames(mMarginal)
str(dfData)
rownames(dfCluster.name) = dfCluster.name$V2
dfCluster.name[cn,]
## create short sensible names for the clusters manually to plot the data
## look at the saved mapping of cluster to gene list table 


cn = c('Generic Transcription',
       'Metabolism of lipids',
       'Transport of small mol',
       'Vesicle-mediated transport',
       'Metabolism of amino acids',
       'Innate Immune System',
       'Pre-mRNA Processing',
       'Biological oxidations',
       'Ext-cell matrix organization',
       'Cell-Cell communication',
       'Hemostasis',
       'Axon guidance',
       'Metabolism of RNA',
       'Adaptive Immune System',
       'Phospholipid metabolism',
       'Signaling by GPCR',
       'Signal Receptor Tyr Kinase',
       'Neuronal System',
       'rRNA processing',
       'Keratinization',
       'White adipocyte diff',
       'Carbohydrate metabolism',
       'Response external stimuli',
       'Deubiquitination',
       'Vitamins, cofactors metab',
       'Cytokine Signaling',
       'GPI-anchored proteins',
       'Cell cycle',
       'Translation'
       )
colnames(dfData) = cn

## sanity check for names
temp = dfCluster.name[colnames(mMarginal),]
temp$cn = cn
temp
#colnames(dfData) = gsub('X', '', colnames(dfData))
dfData$groups = factor(rownames(dfData))
dfData$time = factor(gsub('^Control:|Keloid:', '', as.character(dfData$groups)))
dfData$condition = factor(gsub(':1|:2', '', as.character(dfData$groups)))

dfStack = stack(dfData)
dfStack$time = dfData$time
dfStack$condition = dfData$condition

xyplot(values ~ time | ind, data=dfStack, type=c('b'), groups=condition,  par.strip.text=list(cex=0.7),
       scales=list(cex=0.6), xlab='Time', ylab='Scaled Module Average')#auto.key = list(columns=2), layout=c(4,7))

## add cluster names and gene names to the table
i = match(dfCluster$cluster, dfCluster.name$V2)
dfCluster$clusterName = dfCluster.name$V4[i]
df = f_dfGetGeneAnnotation(as.character(dfCluster$gene))
dfCluster = cbind(dfCluster, SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
write.csv(dfCluster, file='results/dfCluster_members.xls')