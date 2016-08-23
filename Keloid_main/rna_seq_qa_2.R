# File: rna_seq_qa_2.R
# Auth: uhkniazi
# DESC: quality checks on the rna-seq data
# Date: 23/08/2016


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
setwd('Data_external/keloid1')
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
fLane = strsplit(csFiles, '_')
fLane = sapply(fLane, function(x) x[3])

# split by samples
fSamples = strsplit(csFiles, '_')
fSamples = sapply(fSamples, function(x) x[1])

table(fLane, fSamples)

# add sample names and file location
f = fSamples %in% dfSample$title
table(f)

i = match(fSamples, dfSample$title)
dfFiles = data.frame(fSamples, title=dfSample$title[i], files=csFiles, fLane)

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

n = 1:nrow(dfFiles)
lOb = lapply(n, function(x) {
  write.qa(as.character(dfFiles[x,'files']), indir=getwd(), title=paste(dfFiles[x,'title'], dfFiles[x, 'fLane']))
})

setwd(gcswd)
dir.create('Keloid_main/Objects')
save(lOb, file='Keloid_main/Objects/CFastqQuality.object.s021.rds')

getwd()
pdf(file='Keloid_main/Results/qa.quality.cycle.s021.pdf')
par(mfrow=c(2,1))
temp = sapply(lOb, plot.qualitycycle)
dev.off(dev.cur())

