# File: gsea_summary_tables.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 02/12/2016
# Desc: summary tables after merging gsea data

## set variables and source libraries
gcswd = getwd()
setwd('Keloid_main/')
p.old = par()

lFiles = list.files('Results/', pattern='*pathways_mSigDb*', full.names = T, ignore.case = F)

# load the files
ldfData = lapply(lFiles, function(x) read.csv(x, row.names=1))
names(ldfData) = lFiles
sapply(ldfData, nrow)
# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

# get the upregulated/downregulated in order
ldfData.up = ldfData[grepl(lFiles, pattern = 'upregulated')]
ldfData.down = ldfData[grepl(lFiles, pattern = 'downregulated')]

## set the names for each contrast
sn = gsub('Results//(\\w+:\\w+:\\w)_\\w+\\.xls', '\\1', names(ldfData.up))
names(ldfData.up) = sn

sn = gsub('Results//(\\w+:\\w+:\\w)_\\w+\\.xls', '\\1', names(ldfData.down))
names(ldfData.down) = sn

## create a table/matrix of p-values
mMerged.up = sapply(ldfData.up, function(x) x$p.val)
rownames(mMerged.up) = rownames(ldfData.up[[1]])
## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])

# sanity check
identical(rownames(mMerged.up), rownames(mMerged.down))
cn = c('C2vsC1', 'K1vsC1', 'K2vsC2', 'K2vsK1')
colnames(mMerged.up) = paste(cn, 'up', sep='-')
colnames(mMerged.down) = paste(cn, 'down', sep='-')

mMerged.c2 = cbind(mMerged.up, mMerged.down)
# reorder the columns
colnames(mMerged.c2)
o = c(1, 5, 2, 6, 3, 7, 4, 8)
mMerged.c2 = mMerged.c2[,o]

# remove na sections
dim(mMerged.c2)
mMerged.c2 = na.omit(mMerged.c2)
dim(mMerged.c2)
head(mMerged.c2)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) round(x, 3) <= cutoff)
}

mMerged.c2.bin = getBinaryMatrix(mMerged.c2)

## group this matrix into combinations
mMerged.c2.bin.grp = mMerged.c2.bin
set.seed(123)
dm = dist(mMerged.c2.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.c2.bin.grp = cbind(mMerged.c2.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.c2.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## map these names to the cp
groups = cp
sig.pvals = rowSums(mMerged.c2.bin)
dfMerged.c2 = data.frame(round(mMerged.c2, 3), sig.pvals, groups, DB='C2')
str(dfMerged.c2)
head(dfMerged.c2)
tail(dfMerged.c2)

### repeat for the C5 database
lFiles = list.files('Results/', pattern='*msigdb_c5*', full.names = T, ignore.case = F)

# load the files
ldfData = lapply(lFiles, function(x) read.csv(x, row.names=1))
names(ldfData) = lFiles
sapply(ldfData, nrow)
# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

# get the upregulated/downregulated in order
ldfData.up = ldfData[grepl(lFiles, pattern = 'upregulated')]
ldfData.down = ldfData[grepl(lFiles, pattern = 'downregulated')]

## set the names for each contrast
sn = c('C2vsC1-up', "K1vsC1-up", "K2vsK1-up", "K2vsC2-up")
names(ldfData.up) = sn

sn = c("C2vsC1-down", "K1vsC1-down",  "K2vsK1-down", "K2vsC2-down")
names(ldfData.down) = sn

## create a table/matrix of p-values
mMerged.up = sapply(ldfData.up, function(x) x$p.val)
rownames(mMerged.up) = rownames(ldfData.up[[1]])
## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])

# sanity check
identical(rownames(mMerged.up), rownames(mMerged.down))

mMerged.c5 = cbind(mMerged.up, mMerged.down)
# reorder the columns
colnames(mMerged.c2)
colnames(mMerged.c5)
o = match(colnames(mMerged.c2), colnames(mMerged.c5))
mMerged.c5 = mMerged.c5[,o]
identical(colnames(mMerged.c2), colnames(mMerged.c5))

# remove na sections
dim(mMerged.c5)
mMerged.c5 = na.omit(mMerged.c5)
dim(mMerged.c5)
head(mMerged.c5)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) round(x, 3) <= cutoff)
}

mMerged.c5.bin = getBinaryMatrix(mMerged.c5)

## group this matrix into combinations
mMerged.c5.bin.grp = mMerged.c5.bin
set.seed(123)
dm = dist(mMerged.c5.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.c5.bin.grp = cbind(mMerged.c5.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.c5.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## map these names to the cp
groups = cp
sig.pvals = rowSums(mMerged.c5.bin)
table(sig.pvals)
dfMerged.c5 = data.frame(round(mMerged.c5, 3), sig.pvals, groups, DB='C5')
str(dfMerged.c5)
head(dfMerged.c5)
tail(dfMerged.c5)

identical(colnames(dfMerged.c2), colnames(dfMerged.c5))

# #### repeat for the downregulated binary groupings
# mMerged.down.bin.grp = mMerged.down.bin
# set.seed(123)
# dm = dist(mMerged.down.bin.grp, method='binary')
# hc = hclust(dm)
# plot(hc, labels=F)
# 
# # cut the tree at the bottom to create groups
# cp = cutree(hc, h = 0.2)
# # sanity checks
# table(cp)
# length(cp)
# length(unique(cp))
# 
# mMerged.down.bin.grp = cbind(mMerged.down.bin.grp, cp)
# 
# ### print and observe this table and select the groups you are interested in
# temp = mMerged.down.bin.grp
# temp = (temp[!duplicated(cp),])
# temp2 = cbind(temp, table(cp))
# rownames(temp2) = NULL
# print(temp2)
# 
# # create names for the groups
# cvGroupNames = c('C2vsC1-K1vsC1', 'C2vsC1', 'NONE', 'K1vsC1', 'K2vsK1', 'K2vsC2', 'K1vsC1-K2vsC2')
# rownames(temp2) = cvGroupNames
# 
# ## map these names to the cp
# i = match(cp, temp2[,'cp'])
# groups = rownames(temp2)[i]
# sig.pvals = rowSums(mMerged.down.bin)
# dfMerged.down = data.frame(round(mMerged.down, 3), sig.pvals, groups)
# str(dfMerged.down)
# head(dfMerged.down)
# tail(dfMerged.down)
# 
# lVenn = lapply(colnames(mMerged.down.bin), function(n){
#   x = mMerged.down.bin[,n]
#   return(rownames(mMerged.down.bin)[x])
# })
# 
# names(lVenn) = c('C2vsC1', 'K1vsC1', 'K2vsC2', 'K2vsK1')
# venn.diagram(lVenn, filename = 'Results/venn_all_downregulated_msigdb.tif')
# 
write.csv(dfMerged.c2, file='Results/gsea_msigdb_c2_merged.xls')
write.csv(dfMerged.c5, file='Results/gsea_msigdb_c5_merged.xls')

## merge together into one dataframe
table(dfMerged.c2$groups)
dfMerged.c2.sub = dfMerged.c2[dfMerged.c2$groups != 3,]

table(dfMerged.c5$groups)
dfMerged.c5.sub = dfMerged.c5[dfMerged.c5$groups != 2,]

dfMerged = rbind(dfMerged.c2.sub, dfMerged.c5.sub)
dfMerged = droplevels.data.frame(dfMerged)

write.csv(dfMerged, file='Results/gsea_msigdb_significant_c2_c5_merged.csv')

### heatmaps
df = dfMerged
head(df)
mMat = as.matrix(df[,c(1:8)])
head(mMat)
mMat = -10*log10(mMat+1e-16)
g1 = df[,'groups']
g1 = factor(as.character(g1))
levels(g1)
g2 = df[,'DB']
g2 = factor(as.character(g2))
levels(g2)

ann = data.frame(DB=g2, Group=g1 )
range(mMat)
mMat[mMat < 15] = 0 
mMat[mMat > 100] = 100

library(NMF)
library(RColorBrewer)
aheatmap(mMat, annRow = ann, scale = 'none', Rowv = order(g2:g1), Colv=NA, cexRow=5, cexCol = 0.6, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))

pdf('Results/gsea_msigdb_significant_merged.pdf')
aheatmap(mMat, annRow = ann, scale = 'none', Rowv = order(g2:g1), Colv=NA, cexRow=5, cexCol = 0.6, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))
dev.off(dev.cur())
