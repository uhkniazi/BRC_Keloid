# File: counts_from_bams.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: generate count tables for transcripts from bam files
# Date: 15/09/2016


## set variables and source libraries
gcswd = getwd()
setwd('Keloid_main/')

## load the transcript db objects
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicAlignments)

# get the exons into GRangesList object
oGRLgenes = exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = 'gene')

## create the bamfile list from database
## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'File')
# get the query
q = paste0('select Sample.id as sid, Sample.group1 as timepoints, Sample.title, File.* from Sample, File
where (Sample.idData = 2) AND (Sample.group1 like "%Timepoint%") AND (File.idSample = Sample.id AND File.type like "%duplicates removed%")')
dfSample = dbGetQuery(db, q)
# close connection after getting data
dbDisconnect(db)
# remove any whitespace from the names
dfSample$name = gsub(" ", "", dfSample$name, fixed = T)
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)
dfSample$group1 = gsub(" ", "", dfSample$group1, fixed = T)
# create a new file path variable as each sequencing run is in a different directory
dfSample$fp = paste0(dfSample$group1, '/', dfSample$name)

#### set working directory to appropriate location with bam files
setwd('Data_external/Aligned/')
csFiles = list.files('.', pattern = '*.bam', recursive = T)
# check if these files match the file names in database
table(dfSample$fp %in% csFiles)

csFiles = dfSample$fp

## create a bamfiles list object
oBamFiles = BamFileList(csFiles)

oSummExp = summarizeOverlaps(oGRLgenes, oBamFiles, ignore.strand = F, singleEnd=F)

## save the summarized experiment object
setwd(gcswd)
n = make.names(paste('oSummarized experiment counts object for keloids rds'))
n2 = paste0('~/Data/MetaData/', n)
save(oSummExp, file=n2)

library('RMySQL')
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
df = data.frame(idData=2, name=n, type='rds', location='~/Data/MetaData/', 
                comment='oSummarized experiment counts object for keloids S014 S021 and S032 sequencing runs with quality 10 duplicates removed')
dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)

