# File: 01_createDbEntry.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: list the samples and create appropriate db entries
# Date: 21/1/2020


## set variables and source libraries
source('header.R')

## connect to mysql database 
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# sample and file table
dbGetQuery(db, paste('describe Sample;'))
cSampleCol = dbGetQuery(db, paste('describe Sample;'))$Field[-1]

dbGetQuery(db, paste('describe File;'))
cFileCol = dbGetQuery(db, paste('describe File;'))$Field[-1]

# setwd(gcRemoteDir)
setwd('dataExternal/')
setwd('remoteData/raw/')

# list the files
cvFiles = list.files(pattern = 'fastq.gz', recursive = T)

# data is single end 
# # each sample has 2 files 
# fSplit = gsub('_R[1|2]', '', cvFiles)
# 
# lFiles = split(cvFiles, fSplit)

setwd(gcswd)
## load the metadata file
dfMeta = read.csv('dataExternal/metaData.csv', header=T, stringsAsFactors = F)
str(dfMeta)

# sanity check
table(as.character(dfMeta$File_name) %in% cvFiles)

## order the table in the same sequence as file names
i = match(cvFiles, dfMeta$File_name)
dfMeta = dfMeta[i,]
identical(as.character(dfMeta$File_name), cvFiles)
#dfMeta$fSplit = fSplit

## extract reduced sample table by removing one pair of the fastq file
# cvSample = unique(dfMeta$fSplit)
# i = match(cvSample, dfMeta$fSplit)
# dfMeta.sam = dfMeta[i,]
dfMeta.sam = dfMeta
str(dfMeta.sam)
xtabs( ~ Genotype + Litter_ID, data=dfMeta.sam)

## create the entry for samples
cSampleCol

dfSamples = data.frame(idProject=g_pid, idData=g_did, title=dfMeta.sam$File_name, 
                       description= paste('group1 is Treatment or Genotype',
                                          'group2 is Litter_ID',
                                          'group3 is sequencing lane',
                                          sep=';'),
                       group1 = dfMeta.sam$Genotype, group2= dfMeta.sam$Litter_ID, group3=dfMeta.sam$Lane)
# write this data to the database
rownames(dfSamples) = NULL

### NOTE: Do not execute this anymore as entry created
# write this table to database
#dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)
# get this table again from database with ids added
g_did
dfSamples = dbGetQuery(db, paste('select * from Sample where Sample.idData = 44;'))

# create entries for these files in the database
dbListTables(db)
cn = dbListFields(db, 'File')[-1]
cn
# get the names of the samples
# temp = lapply(dfSamples$title, function(x){
#   # get the file names
#   df = data.frame(name=lFiles[[x]], type='fastq', idSample=dfSamples[dfSamples$title == x, 'id'])
#   return(df)
# })

#dfFiles = do.call(rbind, temp)
dfFiles = data.frame(name=cvFiles, type='fastq', idSample=dfSamples$id)
identical(cvFiles, dfSamples$title)
rownames(dfFiles) = NULL

# write this table to database
## note: do not execute as it is already done
# dbWriteTable(db, name='File', value=dfFiles, append=T, row.names=F)

dbDisconnect(db)
