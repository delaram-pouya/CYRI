
Initialize = function()
{
  #### Loading required libraries ####
  options(stringsAsFactors = F)
  set.seed(123) 
  
  #BiocManager::install('pacman')
  
  pac <- list( 'ggplot2', 'limma', 'pheatmap','Biostrings','dplyr','GenomicRanges',
               'VennDiagram','e1071' ,'reshape2', 'ggrepel', 'RColorBrewer', 'caret','ggfortify',
               'plyr', 'gridExtra','caTools', 'h2o','gtools','stringr', "pROC", 'ROCR', 'FactoMineR',
               'factoextra', 'data.table', 'wesanderson', 'viridis')
  
  print(paste(pac , lapply(pac, require, character.only = TRUE), sep = ' : '))
  
  pacman::p_load( 'ggplot2', 'limma', 'pheatmap', 'Biostrings','dplyr','GenomicRanges',
                  'VennDiagram','e1071' ,'reshape2', 'ggrepel', 'RColorBrewer', 'caret','ggfortify',
                  'plyr','gridExtra','caTools','h2o','gtools','stringr', "pROC", 'ROCR', 'FactoMineR',
                  'factoextra','data.table', 'wesanderson', 'viridis')
}

## functions
.cleanHub <- function(x){
  x <- x[-1, ]
  names(x) <- x[1, ]
  x <- x[-1, ]
  colnames(x)[1:3] <- c('index', 'id', 'Gene')
  return(x)
}

.sortPair <- function(i, DF){
  aPair <- sort(DF[i, c('ADORF' ,'DBORF')])
  paste0(aPair[1],'_' ,aPair[2])
}

.getPair <- function(x) sapply(1:nrow(x), function(i) .sortPair(i, x))

.getAD_Degrees <- function(x) {
  AD_deg <- data.frame(table(x$ADORF))
  AD_deg <- AD_deg[order(AD_deg$Freq, decreasing = T),]
  colnames(AD_deg) <- c('Gene', 'Freq'); AD_deg
}

.getDB_Degrees <- function(x) {
  DB_deg <- data.frame(table(x$DBORF))
  DB_deg <- DB_deg[order(DB_deg$Freq, decreasing = T),]
  colnames(DB_deg) <- c('Gene', 'Freq'); DB_deg
}
