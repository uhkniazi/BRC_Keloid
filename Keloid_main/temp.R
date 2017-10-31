dfContrast1 = read.csv(file='Results/Control:2vsControl:1.csv')
dfContrast2 = read.csv(file='Results/Keloid:1vsControl:1.csv')
dfContrast3 = read.csv(file='Results/Keloid:2vsKeloid:1.csv')
dfContrast4 = read.csv(file='Results/Keloid:2vsControl:2.csv')

getFalsePositives = function(x, df){
   nr = nrow(na.omit(df[df$p.adj <= x & df$Dispersion > 0.4,]))
   fp = nr*x
   return(fp)
}

ldf = list(dfContrast1, dfContrast2, dfContrast3, dfContrast4)
cutoffs = seq(0.01, 1, by = 0.02)
mCuts = sapply(ldf, function(x){
  yval = sapply(cutoffs, getFalsePositives, x)
  return(yval)
})

rownames(mCuts) = cutoffs
colnames(mCuts) = c('C1', 'C2', 'C3', 'C4')

matplot(mCuts, type='l', lty=1, xlab='FDR Cutoffs', ylab='False Positives', main='No. of false positives and FDR threshold', col=1:4, lwd=2,
        xaxt='n')
axis(1, at = 1:length(cutoffs), labels = cutoffs, las=2, cex.axis=0.6)
legend('topleft', legend = colnames(mCuts), fill=1:4)

# curve(getFalsePositives, 0.01, 0.2,n = 10, type = 'l', xlab = 'P.adjust cutoff', ylab='False positives', main='Contrast 1')
# 
# df = dfContrast2
# curve(getFalsePositives, 0.01, 1,n = 100, type = 'l', xlab = 'P.adjust cutoff', ylab='False positives', main='Contrast 2')

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
oGRLgenes = exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = 'gene')
i = i[names(i) %in% names(oGRLgenes)]
oGRLgenes = oGRLgenes[names(i)]
f = sapply(width(oGRLgenes), sum) > 3000
oGRLgenes = oGRLgenes[f]
length(oGRLgenes)

# ### take a sample of exons/transcripts
# oGRLgenes = exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = 'tx')
# ## select transcripts with length close to  mean
oGRLgenes = oGRLgenes[sample(1:length(oGRLgenes), 3000, replace = F)]
length(oGRLgenes)
oGRLexons = oGRLgenes
# temp = unlist(oGRLgenes)
# temp = temp[width(temp) > 120]
# temp = temp[sample(1:length(temp), size = 10000, replace = F)]
# oGRLexons = split(temp, 1:length(temp))

getTranscriptCoverage = function(bam){
  ###
  f_bin_vector = function(start, end, bins){
    s = floor(seq(start, end, length.out=bins+1))
    e = s-1
    e[length(e)] = s[length(s)]
    length(s) = length(s)-1
    e = e[2:length(e)]
    return(data.frame(start=s, end=e))
  }# f_bin_vector
  ###
  which = unlist(range(oGRLgenes))
  which = resize(which, width = width(which)+100, fix='center')
  param = ScanBamParam(flag=scanBamFlag(), what = scanBamWhat(), which=which)
  # read the GAlignments object
  oGA = readGAlignments(bam, param=param)
  # get the coverage
  cov = coverageByTranscript(oGA, oGRLexons, ignore.strand=FALSE)
  mCoverage = sapply(cov, function(temp) {
    # create a binned vector to create views on this coverage
    bins = f_bin_vector(1, sum(width(temp)), bins=2000)
    # create views on these bins
    ivCoverage = viewMeans(Views(temp, bins$start, bins$end))
    ivCoverage = ivCoverage / max(ivCoverage)})
  mCoverage = t(mCoverage)
  mCoverage = colMeans(mCoverage, na.rm = T)
  mCoverage = mCoverage/max(mCoverage)
  return(mCoverage)
}
