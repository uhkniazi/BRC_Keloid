# File: gsea_summary_tables.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 02/12/2016
# Desc: summary tables after merging gsea data

## set variables and source libraries
gcswd = getwd()
setwd('Keloid_main/')
p.old = par()

lFiles = list.files('Results/', pattern='*pathways_mSigDb*', full.names = T)

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
# remove na sections
dim(mMerged.up)
mMerged.up = na.omit(mMerged.up)
dim(mMerged.up)
head(mMerged.up)


## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])
# remove na sections
dim(mMerged.down)
mMerged.down = na.omit(mMerged.down)
dim(mMerged.down)
head(mMerged.down)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) x < cutoff)
}

mMerged.up.bin = getBinaryMatrix(mMerged.up)
mMerged.down.bin = getBinaryMatrix(mMerged.down)

## group this matrix into combinations
mMerged.up.bin.grp = mMerged.up.bin
set.seed(123)
dm = dist(mMerged.up.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.up.bin.grp = cbind(mMerged.up.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.up.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

# create names for the groups
cvGroupNames = c('NONE', 'K2vsC2', 'K1vsC1', 'K1vsC1-K2vsC2', 'K2vsK1', 'K2vsC2-K2vsK1', 'C2vsC1', 'C2vsC1-K1vsC1', 
                 'C2vsC1-K2vsK1', 'C2vsC1-K2vsC2-K2vsK1', 'ALL', 'C2vsC1-K1vsC1-K2vsK1')
rownames(temp2) = cvGroupNames

## map these names to the cp
i = match(cp, temp2[,'cp'])
groups = rownames(temp2)[i]
sig.pvals = rowSums(mMerged.up.bin)
dfMerged.up = data.frame(round(mMerged.up, 3), sig.pvals, groups)
str(dfMerged.up)
head(dfMerged.up)
tail(dfMerged.up)

library(VennDiagram)
lVenn = lapply(colnames(mMerged.up.bin), function(n){
  x = mMerged.up.bin[,n]
  return(rownames(mMerged.up.bin)[x])
})

names(lVenn) = c('C2vsC1', 'K1vsC1', 'K2vsC2', 'K2vsK1')
venn.diagram(lVenn, filename = 'Results/venn_all_upregulated_msigdb.tif')

#### repeat for the downregulated binary groupings
mMerged.down.bin.grp = mMerged.down.bin
set.seed(123)
dm = dist(mMerged.down.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)

# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.down.bin.grp = cbind(mMerged.down.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.down.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

# create names for the groups
cvGroupNames = c('C2vsC1-K1vsC1', 'C2vsC1', 'NONE', 'K1vsC1', 'K2vsK1', 'K2vsC2', 'K1vsC1-K2vsC2')
rownames(temp2) = cvGroupNames

## map these names to the cp
i = match(cp, temp2[,'cp'])
groups = rownames(temp2)[i]
sig.pvals = rowSums(mMerged.down.bin)
dfMerged.down = data.frame(round(mMerged.down, 3), sig.pvals, groups)
str(dfMerged.down)
head(dfMerged.down)
tail(dfMerged.down)

lVenn = lapply(colnames(mMerged.down.bin), function(n){
  x = mMerged.down.bin[,n]
  return(rownames(mMerged.down.bin)[x])
})

names(lVenn) = c('C2vsC1', 'K1vsC1', 'K2vsC2', 'K2vsK1')
venn.diagram(lVenn, filename = 'Results/venn_all_downregulated_msigdb.tif')

write.csv(dfMerged.up, file='Results/gsea_msigdb_upregulated_merged.xls')
write.csv(dfMerged.down, file='Results/gsea_msigdb_downregulated_merged.xls')

