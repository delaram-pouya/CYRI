## in this file we'll try to clean the annotation file
source('Code/Functions.R')
Initialize()

annotationFile <- read.delim('Data/Sources/genes.gtf', header = F)
annotSeq <- read.delim('Data/Sources/genomeNew.fa')

split_col <- str_split(annotationFile$V9, '; ')
table(unlist(lapply(split_col, length)))  ## some rows contain TSS-id and ? -> we'll ignore those columns
temp_clean_df <- data.frame(sapply(1:6, function(i) unlist(lapply(split_col, function(x) x[i]))))
column_names <- unlist(lapply(temp_clean_df, function(x) str_split_fixed(x, ' ', 2)[1]))
listOfCols <- lapply(temp_clean_df, function(x) str_split_fixed(x, ' ', 2)[,2])
extraInfo <- data.frame(do.call(cbind, listOfCols))
colnames(extraInfo) <- column_names
head(extraInfo)
annotationFile_clr <- cbind(annotationFile[,1:8], extraInfo)
saveRDS(annotationFile_clr, 'Data/Cleaned_annotation_file.rds')

