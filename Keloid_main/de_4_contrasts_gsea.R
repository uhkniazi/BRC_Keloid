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

