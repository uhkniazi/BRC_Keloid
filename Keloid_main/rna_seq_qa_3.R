# File: rna_seq_qa_3.R
# Auth: uhkniazi
# DESC: quality checks on the rna-seq data
# Date: 08/09/2016


## set variables and source libraries
gcswd = getwd()
setwd('Keloid_main/')
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CFastqQuality/master/CFastqQuality.R'
download(url, 'CFastqQuality.R')

# load the required packages
source('CFastqQuality.R')
# delete the file after source
unlink('CFastqQuality.R')

## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'Sample')
# another way to get the query, preferred
dfSample = dbGetQuery(db, "select title, group1, group2 from Sample where idData=2;")
# close connection after getting data
dbDisconnect(db)
# remove any whitespace from the names
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)

#### get the names of the fastq files for first sequencing run
setwd('Data_external/Trimmed/S032/')
csFiles = list.files('.', pattern = '*.gz')

# match the files to samples by title
i = sapply(dfSample$title, function(x) grep(x, csFiles))

# index of list with no matching files
i2 = sapply(i, length)
# skip this step if i2 has no 0 values
dfSample = dfSample[which(i2 != 0),]

# match the files to their respective index in the table
lFilesIndex = lapply(dfSample$title, function(x) grep(x, csFiles))
names(lFilesIndex) = dfSample$title


## perform the analysis one sample at a time
## function to write the qa files
write.qa = function(fls, indir, title){
  wd = getwd()
  setwd(indir)
  ob = CFastqQuality(fls, title)
  setwd(wd)
  cat(paste('done', title, '\n'))
  return(ob)
}

ivFilesIndex = unlist(lFilesIndex)

lOb = lapply(ivFilesIndex, function(x){
  write.qa(csFiles[x], getwd(), csFiles[x])
})

setwd(gcswd)
n = make.names(paste('CFastqQuality trimmed S032 rds'))
lOb$meta = dfSample
lOb$desc = paste('CFastqQuality object from S032 sequencing run', date())
n2 = paste0('~/Data/MetaData/', n)
save(lOb, file=n2)

library('RMySQL')
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
df = data.frame(idData=2, name=n, type='rds', location='~/Data/MetaData/', comment='Sequencing run S032 FASTQ file quality data after trimming')
dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)

### create the plots of interest
getwd()
lOb$desc = NULL
lOb$meta = NULL
pdf(file='Keloid_main/Results/qa.trimmed.S032.pdf')

iReadCount = sapply(lOb, CFastqQuality.getReadCount)
iReadCount = iReadCount/1e+6

barplot(iReadCount, las=2, main='Sequencing run S032 - Read Count', ylab = 'No. of Reads in Millions', cex.names =0.8, col=grey.colors(2))

mQuality = sapply(lOb, function(x){
  m = mGetReadQualityByCycle(x)
  m = colMeans(m, na.rm = T)
  return(m)
})

matplot(mQuality, type='l', main='Sequencing run S032 - Per base quality', ylab = 'Mean Score', xlab='Position in Read')

lReadWidth = lapply(lOb, iGetReadWidth)
boxplot(lReadWidth, las=2, main='Sequencing run S032 - Read Width', ylab = 'Read Width', cex.names =0.8, col=grey.colors(2))

dev.off(dev.cur())

