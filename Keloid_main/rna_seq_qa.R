# File: rna_seq_qa.R
# Auth: uhkniazi
# DESC: quality checks on the rna-seq data
# Date: 18/08/2016


## source the quality check file
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

#### get the names of the fastq files
setwd('Data_external/keloid1')
csFiles = list.files('.', pattern = '*.gz')

# path to fastq file
csPath = file.choose()

ob = CFastqQuality(csPath, sample.name = 'SRR850132_2.fastq')

plot.alphabetcycle(ob)
plot.qualitycycle(ob)
plot.recurrentreads(ob)

l = lBlastRecurrentSequences(ob, n=2, timeout=100000)


fls <- dir(path = '../temp/RNASeq/20151110/FASTQ/', pattern = "*fastq.gz", full=TRUE)

coll <- QACollate(QAFastqSource(fls), QAReadQuality(),
                  QAAdapterContamination(), QANucleotideUse(),
                  QAQualityUse(), QASequenceUse(),
                  QAFrequentSequence(n=10), QANucleotideByCycle(),
                  QAQualityByCycle())
x <- qa2(coll,  verbose=TRUE)

res <- report(x, dest = 'Temp')
if (interactive())
  browseURL(res)