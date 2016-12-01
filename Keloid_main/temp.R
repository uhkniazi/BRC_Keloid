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