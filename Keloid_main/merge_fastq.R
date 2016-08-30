# File: merge_fastq.R
# Auth: uhkniazi
# DESC: merge fastq files for each sample
# Date: 26/08/2016


## set variables and source libraries
gcswd = getwd()
setwd('Keloid_main/')

## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# check how many files each sample has
q = paste0('select count(File.idSample) as files, Sample.idData, Sample.title, Sample.id from File, Sample
           where (Sample.idData = 2 and File.idSample = Sample.id) group by File.idSample')
dfQuery = dbGetQuery(db, q)
dfQuery$title = gsub(" ", "", dfQuery$title, fixed = T)

# for each sample id, get the corresponding files
cvQueries = paste0('select File.*, Sample.title from File, Sample 
           where (Sample.idData = 2 and Sample.id =', dfQuery$id, ') and (File.idSample = Sample.id)')

#### get the names of the fastq files for S021 or S014 sequencing run
setwd('Data_external/keloid2')
csFiles = list.files('.', pattern = '*.gz')

temp = sapply(cvQueries, function(x){
  # get the file names
  dfFiles = dbGetQuery(db, x)
  # remove white space from title
  dfFiles$title = gsub(" ", "", dfFiles$title, fixed=T)
  f = csFiles[csFiles %in% dfFiles$name]
  # split the file names into paired end 1 and 2, identified by R1 and R2 in the file name
  d = grepl('_R1_', f)
  d = as.character(d)
  d[d == 'TRUE'] = 'R1'
  d[d == 'FALSE'] = 'R2'
  lf = split(f, d)
  # merge the files from each direction
  sapply(names(lf), function(x2){
    # collapse from string vector to character string
    f2 = paste(lf[[x2]], collapse = ' ')
    f.out = paste0(unique(dfFiles$title), '_S014_', x2, '_', '.fastq.gz')
    com = paste('cat', f2, '>>' , f.out)
    cat(com, '\n')
    system(com, intern = F, ignore.stdout = F)
  })
})

# create new entries in the File table for these new files associated with samples
dbListTables(db)
cn = dbListFields(db, 'File')[-1]
# get the names of the samples
temp = lapply(cvQueries, function(x){
  # get the file names
  dfFiles = dbGetQuery(db, x)
  # remove white space from title
  dfFiles$title = gsub(" ", "", dfFiles$title, fixed=T)
  f = csFiles[csFiles %in% dfFiles$name]
  # split the file names into paired end 1 and 2, identified by R1 and R2 in the file name
  d = grepl('_R1_', f)
  d = as.character(d)
  d[d == 'TRUE'] = 'R1'
  d[d == 'FALSE'] = 'R2'
  lf = split(f, d)
  # merge the files from each direction
  name = sapply(names(lf), function(x2){
    f.out = paste0(unique(dfFiles$title), '_S014_', x2, '_', '.fastq.gz')
  })
  df = data.frame(idSample=unique(dfFiles$idSample), name, type='fastq', group1='S014')
  return(df)
})

dfNewData = do.call(rbind, temp)
rownames(dfNewData) = NULL

# write this table to database
dbWriteTable(db, name='File', value=dfNewData, append=T, row.names=F)
dbDisconnect(db)

