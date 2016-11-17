# File: de_4_contrasts_gsea.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 16/11/2016
# Desc: gene set enrichment analysis for the datasets

## set variables and source libraries
gcswd = getwd()
setwd('Keloid_main/')
p.old = par()

## libraries to load
library(gage)

## load gene expression data
dfContrast1 = read.csv(file='Results/Control:2vsControl:1.csv', stringsAsFactors = F)
dfContrast2 = read.csv(file='Results/Keloid:1vsControl:1.csv', stringsAsFactors = F)
dfContrast3 = read.csv(file='Results/Keloid:2vsKeloid:1.csv', stringsAsFactors = F)
dfContrast4 = read.csv(file='Results/Keloid:2vsControl:2.csv', stringsAsFactors = F)

## load msig db data
oMsigGS = readList(file.choose())
## OR try other datasets
library(gageData)
data("kegg.sets.hs")
oMsigGS = kegg.sets.hs
## GO
data(go.sets.hs)
data(go.subs.hs)
names(go.subs.hs)
oMsigGS=go.sets.hs[go.subs.hs$MF]
oMsigGS=go.sets.hs[go.subs.hs$BP]
oMsigGS=go.sets.hs[go.subs.hs$CC]

dfContrast = dfContrast4
# for a contrats of choice create the list
iContFc = dfContrast$logfc
names(iContFc) = as.character(dfContrast[,1])
head(iContFc)
head(dfContrast)

oGage = gage(iContFc, oMsigGS)

dfGreater = data.frame(oGage$greater)
str(dfGreater)
i = which(dfGreater$p.val < 0.01)
rownames(dfGreater[i,])

dfLess = data.frame(oGage$less)
str(dfLess)
i = which(dfLess$p.val < 0.01)
rownames(dfLess[i,])

# write.csv(dfGreater[,c('p.val', 'q.val', 'set.size')], file='Results/Keloid:2vsControl:2_upregulated_pathways_mSigDb_c2_curated.xls')
# write.csv(dfLess[,c('p.val', 'q.val', 'set.size')], file='Results/Keloid:2vsControl:2_downregulated_pathways_mSigDb_c2_curated.xls')

write.csv(dfGreater[,c('p.val', 'q.val', 'set.size')], file='Results/Keloid:2vsControl:2_upregulated_pathways_GO_BP.xls')
write.csv(dfLess[,c('p.val', 'q.val', 'set.size')], file='Results/Keloid:2vsControl:2_downregulated_pathways_GO_BP.xls')



## make a barplot
plot.bar = function(ivBar, title='', ...){
  p.old = par(mar=c(8,4,2,2)+0.1)
  l = barplot(ivBar, beside=T, xaxt='n', main=title, ...)
  axis(side = 1, l[,1], labels=F)
  text(l[,1], y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
       labels=names(ivBar), srt=45, adj=1, xpd=TRUE, cex=0.6)
  par(p.old)
}

dfBar = dfGreater[order(dfGreater$p.val), ]
dfBar = dfLess[order(dfLess$p.val), ]

i = -1*log10(dfBar$p.val[1:7])
names(i) = rownames(dfBar)[1:7]

plot.bar(i, title='Contrast 4 - down regulated GO BP', ylab='-log10 PValue')

#### test with kegg pathways data
library(pathview)
dir.create('Results/Kegg_figures')
pv.out.list <- sapply(substr(rownames(dfBar)[1], 1, 8), function(pid) pathview(
  gene.data = iContFc, pathway.id = pid,
  species = "hsa", out.suffix='Contrast 1 Down', kegg.dir = 'Results/Kegg_figures/'))









### grouping of genes
dfContrast1.sub = na.omit(dfContrast1[dfContrast1$p.adj < 0.01 & dfContrast1$Dispersion > 0.4,])
dfContrast2.sub = na.omit(dfContrast2[dfContrast2$p.adj < 0.1 & dfContrast2$Dispersion > 0.4,])
dfContrast3.sub = na.omit(dfContrast3[dfContrast3$p.adj < 0.01 & dfContrast3$Dispersion > 0.4,])
dfContrast4.sub = na.omit(dfContrast4[dfContrast4$p.adj < 0.1 & dfContrast4$Dispersion > 0.4,])

cvCommonGenes = unique(c(rownames(dfContrast1.sub), rownames(dfContrast2.sub), rownames(dfContrast3.sub), rownames(dfContrast4.sub)))
mCommonGenes = matrix(NA, nrow=length(cvCommonGenes), ncol=4)
mCommonGenes[,1] = cvCommonGenes %in% rownames(dfContrast1.sub)
mCommonGenes[,2] = cvCommonGenes %in% rownames(dfContrast2.sub)
mCommonGenes[,3] = cvCommonGenes %in% rownames(dfContrast3.sub)
mCommonGenes[,4] = cvCommonGenes %in% rownames(dfContrast4.sub)
rownames(mCommonGenes) = cvCommonGenes
colnames(mCommonGenes) = gsub(' ', '', rownames(mContrasts))


############

#### analysis by grouping genes
# create groups in the data based on 4^2-1 combinations
mCommonGenes.grp = mCommonGenes
set.seed(123)
dm = dist(mCommonGenes.grp, method='binary')
hc = hclust(dm)
plot(hc)
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

## make venn diagrams of comparisons
library(VennDiagram)
par(p.old)
# create a list for overlaps
lVenn = list(rownames(dfContrast1.sub), rownames(dfContrast2.sub), rownames(dfContrast3.sub), rownames(dfContrast4.sub))
names(lVenn) = c('C2vsC1', 'K1vsC1', 'K2vsK1', 'K2vsC2')
# calculate overlaps
#lVenn.overlap = calculate.overlap(lVenn)
venn.diagram(lVenn, filename = 'Results/venn_all_contrasts.tif')
venn.diagram(lVenn[c(1,3)], filename = 'Results/venn_time_contrasts.tif')

## save the genes in the overlaps of interest
# genes common between C2vsC1 and K2vsK1
rn = apply(mCommonGenes, 1, function(x) all(x[c(1, 3)] == c(T, T)))
rn = names(rn[rn])
length(rn)

df.rn = select(org.Hs.eg.db, keys = rn, columns = c('SYMBOL', 'GENENAME'), keytype = 'ENTREZID')
dim(df.rn)
str(df.rn)
# write csv to look at gene list
write.csv(df.rn, file=paste('Temp/', 'common.between.c2vsc1.and.k2vsk1', '.xls', sep=''))

# genes unique to C2vsC1 VS K2vsK1
rn = apply(mCommonGenes, 1, function(x) all(x[c(1, 3)] == c(T, F)))
rn = names(rn[rn])
length(rn)

df.rn = select(org.Hs.eg.db, keys = rn, columns = c('SYMBOL', 'GENENAME'), keytype = 'ENTREZID')
dim(df.rn)
str(df.rn)
# write csv to look at gene list
write.csv(df.rn, file=paste('Temp/', 'unique.to.c2vsc1.VS.k2vsk1', '.xls', sep=''))

# genes unique to K2vsK1 VS C2vsC1
rn = apply(mCommonGenes, 1, function(x) all(x[c(1, 3)] == c(F, T)))
rn = names(rn[rn])
length(rn)

df.rn = select(org.Hs.eg.db, keys = rn, columns = c('SYMBOL', 'GENENAME'), keytype = 'ENTREZID')
dim(df.rn)
str(df.rn)
# write csv to look at gene list
write.csv(df.rn, file=paste('Temp/', 'unique.to.k2vsk1.VS.c2vsc1', '.xls', sep=''))

# genes present in other 2 remaining contrasts
df.rn = select(org.Hs.eg.db, keys = rownames(dfContrast2.sub), columns = c('SYMBOL', 'GENENAME'), keytype = 'ENTREZID')
dim(df.rn)
str(df.rn)
write.csv(df.rn, file='Temp/K1vsC1.xls')

df.rn = select(org.Hs.eg.db, keys = rownames(dfContrast4.sub), columns = c('SYMBOL', 'GENENAME'), keytype = 'ENTREZID')
dim(df.rn)
str(df.rn)
write.csv(df.rn, file='Temp/K2vsC2.xls')

## insert this gene list at innate db to see pathways overrepresented
## import innate db result
dfCommonAcrossTime = read.csv(file.choose(), header = T, sep='\t', stringsAsFactors = F)
dfCommonAcrossTime = dfCommonAcrossTime[,-10]
# sort the pathway table on adjusted p-values
dfCommonAcrossTime = dfCommonAcrossTime[order(dfCommonAcrossTime$Pathway.p.value),]

## make a barplot
plot.bar = function(ivBar, title='', ...){
  p.old = par(mar=c(8,4,2,2)+0.1)
  l = barplot(ivBar, beside=T, xaxt='n', main=title, ...)
  axis(side = 1, l[,1], labels=F)
  text(l[,1], y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
       labels=names(ivBar), srt=45, adj=1, xpd=TRUE, cex=0.6)
  par(p.old)
}

i = -1*log10(dfCommonAcrossTime$Pathway.p.value[1:7])
names(i) = dfCommonAcrossTime$Pathway.Name[1:7]
# make on of the names shorter
names(i)[2] = 'Extracellular matrix degradation'
pdf('Temp/injury_response.pdf')
plot.bar(i, title='Common Response to Injury', ylab='-log10 PValue')
dev.off(dev.cur())

## recycling the code above 
## import innate db result for unique to keloids across time
dfCommonAcrossTime = read.csv(file.choose(), header = T, sep='\t', stringsAsFactors = F)
dfCommonAcrossTime = dfCommonAcrossTime[,-10]
# sort the pathway table on adjusted p-values
dfCommonAcrossTime = dfCommonAcrossTime[order(dfCommonAcrossTime$Pathway.p.value),]

## make a barplot
plot.bar = function(ivBar, title='', ...){
  p.old = par(mar=c(8,4,2,2)+0.1)
  l = barplot(ivBar, beside=T, xaxt='n', main=title, ...)
  axis(side = 1, l[,1], labels=F)
  text(l[,1], y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
       labels=names(ivBar), srt=45, adj=1, xpd=TRUE, cex=0.6)
  par(p.old)
}

i = -1*log10(dfCommonAcrossTime$Pathway.p.value[1:7])
names(i) = dfCommonAcrossTime$Pathway.Name[1:7]
# make on of the names shorter
names(i)[7] = 'Molecular Transport'
pdf('Temp/unique_keloids_injury_response.pdf')
plot.bar(i, title='Unique in Keloids Response to Injury', ylab='-log10 PValue')
dev.off(dev.cur())


### make some heatmaps
## all overexpressed genes if interested in
fSamples = oExp$fCondition.t
i = 1:nrow(mCommonGenes)

m1 = as.matrix(mCommonGenes[i,])
m1 = mDat[rownames(m1),]

fGroups = fSamples
colnames(m1) = fGroups
m1 = m1[,order(fGroups)]
fGroups = fGroups[order(fGroups)]

## download the cgraph library
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CGraphClust/master/CGraphClust.R'
download(url, 'CGraphClust.R')

# load the required packages
source('CGraphClust.R')
# delete the file after source
unlink('CGraphClust.R')
library('NMF')
m1 = log(m1+0.5)
# ignore step if stabalization not required
m1 = t(apply(m1, 1, function(x) f_ivStabilizeData(x, fGroups)))
colnames(m1) = fGroups

# scale across rows
dfGenes = select(org.Hs.eg.db, cvCommonGenes, c('SYMBOL', 'GENENAME'), 'ENTREZID')
rownames(dfGenes) = dfGenes$ENTREZID
rownames(m1) = dfGenes[rownames(m1), 'SYMBOL']
mCounts = t(m1)
mCounts = scale(mCounts)
mCounts = t(mCounts)
# threshhold the values
mCounts[mCounts < -3] = -3
mCounts[mCounts > 3] = 3

# draw the heatmap  color='-RdBu:50'
aheatmap(mCounts, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
         annColors=NA, Colv=NA)


####################################################################################################
########### pathway analysis using CGraph library
# uniprot annotation for data
cvGenes = rownames(mCommonGenes)
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
x = gsub('.+\\w+-\\w+-(\\d+)', replacement = '\\1', x = dfReactome$V2, perl = T)
dfReactome$V2 = x
dfReactome.sub = dfReactome[dfReactome$V1 %in% dfGenes$UNIPROT,]
# get the matching positions for uniprot ids in the reactome table
i = match(dfReactome.sub$V1, dfGenes$UNIPROT)
dfReactome.sub$ENTREZID = dfGenes$ENTREZID[i]
dfGraph = dfReactome.sub[,c('ENTREZID', 'V2')]
dfGraph = na.omit(dfGraph)

# get expression data
mCounts = mDat[unique(dfGenes$ENTREZID),]
# ## optional adjusting the data for repeated measurements
# temp = sapply(rownames(mCounts), function(x){
#   return(fitted(lGlm.sub[[x]]))
# })
# mCounts = t(mCounts)

fGroups = fSamples
names(fGroups) = oExp$fTitle
colnames(mCounts) = fGroups
# reorder on grouping factor
mCounts = mCounts[,order(fGroups)]
fGroups = fGroups[order(fGroups)]
mCounts = t(mCounts)

# select genes that have a reactome term attached
n = unique(dfGraph$ENTREZID)
mCounts = mCounts[,n]
print(paste('Total number of genes with Reactome terms', length(n)))
levels(fGroups)

# [1] "Total number of genes with Reactome terms 910"
# > levels(fGroups)
# [1] "Control:1" "Control:2" "Keloid:1"  "Keloid:2" 
# > 

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

write.csv(dfTopGenes.cent, file='Temp/Top_Centrality_Genes.xls')

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
m = getSignificantClusters(oGr, log(t(mCounts)+0.5), fGroups)$clusters
m = getClusterMarginal(oGr, log(t(mCounts)+0.5))
csClust = rownames(m)
length(csClust)
pdf('Temp/cluster_variance.pdf')
par(mfrow=c(1,2))
boxplot.cluster.variance(oGr, m, fGroups, log=T, iDrawCount = length(csClust), las=2)
for (i in seq_along(csClust)){
  temp = t(as.matrix(m[csClust[i],]))
  rownames(temp) = csClust[i]
  plot.cluster.variance(oGr, temp, fGroups, log=FALSE);
}
dev.off(dev.cur())


# print the labels of the clusters from reactome table
i = which(dfReactome.sub$V2 %in% csClust)
dfCluster.name = dfReactome.sub[i,c('V2', 'V4')]
dfCluster.name = dfCluster.name[!duplicated(dfCluster.name$V2),]
rownames(dfCluster.name) = NULL
dfCluster.name

#### plot a graph of clusters 
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
write.csv(dfCluster, file='Temp/dfCluster_members.xls')
# how many genes in each cluster
data.frame(sort(table(dfCluster$cluster)))
#csClust = rownames(m$clusters)
csClust = as.character(unique(dfCluster$cluster))

mMarginal = getClusterMarginal(oGr, log(t(mCounts)+0.5))

### create lattice plots for each cluster and group i.e. keloids and normals/controls
library(lattice)
dfData = data.frame(scale(t(mMarginal)))
#colnames(dfData) = gsub('X', '', colnames(dfData))
dfData$groups = fGroups
dfData$patient = factor(names(fGroups))
# extract patient conditions from the fGroups
fTime = factor(gsub('^Control:|Keloid:', '', as.character(fGroups)))
fCondition = factor(gsub(':1|:2', '', as.character(fGroups)))
dfData$time = fTime
dfData$condition = fCondition

str(dfData)
## stack the dataframe expression values
dfStack = stack(dfData)
str(dfStack)
dfStack$time = dfData$time
dfStack$condition = dfData$condition
dfStack$patient = dfData$patient

xyplot(values ~ time | ind, data=dfStack, type=c('n','smooth'), groups=condition, auto.key=T)

## xyplots with population average
mMarginal = getClusterMarginal(oGr, log(t(mCounts)+0.5))
mMarginal = apply(mMarginal, 1, function(x) tapply((x), fGroups, mean))
mMarginal = scale(mMarginal)
dfData = data.frame(mMarginal)
cn = colnames(mMarginal)
str(dfData)
rownames(dfCluster.name) = dfCluster.name$V2
dfCluster.name[cn,]
cn = c('Transmembrane transport', 'Amino A Met', 'Organelle biogenesis', 'Neutrophil degran/GTPase',
       'Epigenetics', 'ECM Degradation', 'Semaphorin interactions', 'Membrane Transport', 'ECM Proteoglycans',
       'Class I MHC', 'GPCR binding', 'Muscle contraction', 'Sphingolipid metabolism', 'Keratinization',
       'Wnt Signaling', 'Deubiquitination', 'Neuronal System', 'Cell communication', 'Stress response', 'Cytokine Signaling',
       'FA metabolism', 'rRNA processing', 'Lipid metabolism', 'p53 regulation', 'Olf Rec Path')
colnames(dfData) = cn
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
pdf('Temp/module_summary.pdf')
xyplot(values ~ time | ind, data=dfStack, type=c('b'), groups=condition, auto.key=T, par.strip.text=list(cex=0.6),
       scales=list(cex=0.6), xlab='Time', ylab='Scaled Module Average')#, layout=c(4,7))
dev.off(dev.cur())
xyplot(values ~ groups | ind, data=dfStack, type='o')
