# File: rna_seq_qa.R
# Auth: uhkniazi
# DESC: quality checks on the rna-seq data
# Date: 18/08/2016


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
## not the preferred way of calling skip this
# dbSendQuery(db, "select title from Sample where idData=2;")
# # clear last query results
# dbClearResult(dbListResults(db)[[1]])
# this is fine reads whole table
dbReadTable(db, name='Sample')
# another way to get the query, preferred
dfSample = dbGetQuery(db, "select title, group1, group2 from Sample where idData=2;")
# close connection after getting data
dbDisconnect(db)
# remove any whitespace from the names
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)

#### get the names of the fastq files for first sequencing run
setwd('Data_external/keloid2')
csFiles = list.files('.', pattern = '*.gz')

# remove the index files from the list
i = grep('_I', csFiles)
head(csFiles[i])
i2 = file.size(csFiles[i])
head(i2)/1e6
summary(i2/1e6)

csFiles = csFiles[-i]

# remove undetermined
i = grep('Undetermined', csFiles)
head(csFiles[i])
i2 = file.size(csFiles[i])
head(i2)/1e6
summary(i2/1e6)
csFiles = csFiles[-i]

# split samples by lane
# fLane = strsplit(csFiles, '_')
# fLane = sapply(fLane, function(x) x[3])
fLane = gsub('.+(L\\d+)_.+', '\\1', csFiles, perl=T)

# split by samples
# fSamples = strsplit(csFiles, '_')
# fSamples = sapply(fSamples, function(x) x[1])
fSamples = gsub('^([A-Z]+\\d+)_(1st|2nd)_.+', '\\1-\\2', csFiles, perl=T)
fSamples = gsub('^([A-Z]+\\d+)_.+', '\\1', fSamples, perl=T)

table(fLane, fSamples)

# add sample names and file location
f = fSamples %in% dfSample$title
table(f)

i = match(fSamples, dfSample$title)
dfFiles = data.frame(fSamples, title=dfSample$title[i], files=csFiles, fLane)

## perform the analysis one sample at a time
lLanes = split(csFiles, fSamples)
names(lLanes) = unique(fSamples)
sapply(lLanes, length)

## function to write the qa files
write.qa = function(fls, indir, outdir, title){
  wd = getwd()
  setwd(indir)
  coll <- QACollate(QAFastqSource(fls), QAReadQuality(),
                    QAAdapterContamination(), QANucleotideUse(),
                    QAQualityUse(), QASequenceUse(),
                    QAFrequentSequence(n=10), QANucleotideByCycle(),
                    QAQualityByCycle())
  x <- qa2(coll,  verbose=F)
  setwd(outdir)
  res <- report(x, dest = title)
  setwd(wd)
  cat(paste('done', title))
}

n = names(lLanes)
temp = sapply(n, function(x) {
  write.qa(lLanes[[x]], indir=getwd(), outdir='~/Data/R/BRC_Keloid/Keloid_main/Results/', title=x)
})


# write the file information to the database
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'Sample')
## not the preferred way of calling skip this
# another way to get the query, preferred
dfSample = dbGetQuery(db, "select id, title from Sample where idData=2;")
# remove any whitespace from the names
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)
dbGetQuery(db, paste('describe File;'))

# create insert query for each sample
n = names(lLanes)
temp = sapply(n, function(x) {
  # get sample id
  idSample = dbGetQuery(db, paste('select id from Sample where title like "%', x, '%" AND idData=2;', sep=''))
  idSample = cbind(idSample, lLanes[[x]])
  idSample = cbind(idSample, 'fastq')
  i = match(idSample[,2], dfFiles$files)
  idSample$lane = dfFiles$fLane[i]
  colnames(idSample) = c('idSample', 'name', 'type', 'group1')
  ## create a sql query
  dbWriteTable(db, name='File', value=idSample, append=T, row.names=F)
})

# close connection after getting data
dbDisconnect(db)
