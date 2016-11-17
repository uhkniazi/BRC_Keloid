# File: de_4_contrasts_overrepresentation.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 17/11/2016
# Desc: over representation analysis for various overlaps of contrasts

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

rownames(dfContrast1) = dfContrast1$ENTEREZ.ID
rownames(dfContrast2) = dfContrast2$ENTEREZ.ID
rownames(dfContrast3) = dfContrast3$ENTEREZ.ID
rownames(dfContrast4) = dfContrast4$ENTEREZ.ID
### grouping of genes
dfContrast1.sub = na.omit(dfContrast1[dfContrast1$p.adj < 0.01 & dfContrast1$Dispersion > 0.4,])
dfContrast2.sub = na.omit(dfContrast2[dfContrast2$p.adj < 0.1 & dfContrast2$Dispersion > 0.4,])
dfContrast3.sub = na.omit(dfContrast3[dfContrast3$p.adj < 0.01 & dfContrast3$Dispersion > 0.4,])
dfContrast4.sub = na.omit(dfContrast4[dfContrast4$p.adj < 0.1 & dfContrast4$Dispersion > 0.4,])

#### analysis by grouping genes
library(VennDiagram)
par(p.old)
# create a list for overlaps
lVenn = list(as.character(dfContrast1.sub$ENTEREZ.ID), as.character(dfContrast2.sub$ENTREZ.ID),
             as.character(dfContrast3.sub$ENTREZ.ID), as.character(dfContrast4.sub$ENTREZ.ID))
names(lVenn) = c('C2vsC1', 'K1vsC1', 'K2vsK1', 'K2vsC2')

cvCommonGenes = unique(do.call(c, lVenn))
mCommonGenes = matrix(NA, nrow=length(cvCommonGenes), ncol=4)
mCommonGenes[,1] = cvCommonGenes %in% lVenn[['C2vsC1']]
mCommonGenes[,2] = cvCommonGenes %in% lVenn[['K1vsC1']]
mCommonGenes[,3] = cvCommonGenes %in% lVenn[['K2vsK1']]
mCommonGenes[,4] = cvCommonGenes %in% lVenn[['K2vsC2']]
rownames(mCommonGenes) = cvCommonGenes
colnames(mCommonGenes) = names(lVenn)

# create groups in the data based on 4^2-1 combinations
mCommonGenes.grp = mCommonGenes
set.seed(123)
dm = dist(mCommonGenes.grp, method='binary')
hc = hclust(dm)
plot(hc)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.1)
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

## groups of interest
cpi = c(5, 4, 6, 3, 12, 7)

## save the genes in the overlaps of interest
rn = which(mCommonGenes.grp[,'cp'] == cpi[1])
rn = names(rn)
length(rn)

library(org.Hs.eg.db)
df.rn = select(org.Hs.eg.db, keys = rn, columns = c('SYMBOL', 'GENENAME'), keytype = 'ENTREZID')
dim(df.rn)
str(df.rn)
# write csv to look at gene list
write.csv(df.rn, file=paste('Results/', 'venn_overlaps_group_', cpi[1], '.xls', sep=''))

## repeat for other groups
for (i in 2:length(cpi)){
  rn = which(mCommonGenes.grp[,'cp'] == cpi[i])
  rn = names(rn)
  length(rn)
  df.rn = select(org.Hs.eg.db, keys = rn, columns = c('SYMBOL', 'GENENAME'), keytype = 'ENTREZID')
  dim(df.rn)
  str(df.rn)
  # write csv to look at gene list
  write.csv(df.rn, file=paste('Results/', 'venn_overlaps_group_', cpi[i], '.xls', sep=''))
}


## perform innate db over representation analysis and plot results

dfInnat = read.csv(file.choose(), header = T, sep='\t', stringsAsFactors = F)
dfInnat = dfInnat[,-10]
# sort the pathway table on adjusted p-values
dfInnat = dfInnat[order(dfInnat$Pathway.p.value),]

## make a barplot
plot.bar = function(ivBar, title='', ...){
  p.old = par(mar=c(8,4,2,2)+0.1)
  l = barplot(ivBar, beside=T, xaxt='n', main=title, ...)
  axis(side = 1, l[,1], labels=F)
  text(l[,1], y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
       labels=names(ivBar), srt=45, adj=1, xpd=TRUE, cex=0.6)
  par(p.old)
}

i = -1*log10(dfInnat$Pathway.p.value[1:7])
names(i) = dfInnat$Pathway.Name[1:7]
plot.bar(i, title='Venn Overlaps unique to k2 vs k1', ylab='-log10 PValue')

