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

.getComb_Degree <- function(x){
  x <- x[!duplicated(x$PPI), ]  ## removing the PPIs which are duplicated and redundant
  deg <- data.frame(table(c(x$ADORF,x$DBORF)))
  deg <- deg[order(deg$Freq, decreasing = T),]
  colnames(deg) <- c('Gene', 'Freq'); deg
}


getHubComparisonMatrix <- function(ListOfHubs){
  all_hubs <- unlist(lapply(ListOfHubs, function(x) x$Gene))
  print(paste0('#all hubs: ', length(all_hubs)))
  print(paste0('#all unique-hubs: ', length(unique(all_hubs))))
  unique_hubs<- unique(all_hubs)
  compare_hubs<- data.frame(do.call(cbind, lapply(ListOfHubs, function(x) ifelse(unique_hubs %in% x$Gene, 1, 0))))
  compare_hubs$sum <- rowSums(compare_hubs)
  rownames(compare_hubs) <- unique_hubs
  return(compare_hubs)
}


getDegreeMatrix <- function(degrees){
  degrees_tmp <- sapply(1:length(degrees), 
                        function(i) {
                          colnames(degrees[[i]])[2] <- paste0('Freq','_', names(degrees)[i]);degrees[[i]]}, 
                        simplify = F)
  names(degrees_tmp) <- names(degrees)
  degree_Matrix <- Reduce(function(x, y) merge(x, y,by.x='Gene',by.y='Gene' ,all=TRUE), degrees_tmp)
  degree_Matrix[is.na(degree_Matrix)] <- 0
  return(degree_Matrix)
}


get_attributes <- function(input_df, degree_Matrix){
  input_df$id <-rownames(input_df) 
  input_df <- merge(input_df, GOterms, by.x='id' , by.y= 'sys_name', all.x=T)
  merge(input_df, degree_Matrix, by.x='id', by.y= 'Gene', all.x=T, all.y=F)
}



### given a list of tracking IDs , get the corrsponding genes for them
getDEgenes <- function(listOftrackingId){
  DE_genes <- expression[expression$tracking_id %in% listOftrackingId,]$gene_short_name
  DE_genes <- unlist(str_split(DE_genes,','))
  DE_genes <- unique(DE_genes[DE_genes!='-'])
  return(DE_genes)
}


getInteractingPairs <- function(target, ConditionPPIs ){
  c(ConditionPPIs$ADORF[ConditionPPIs$DBORF==target], ConditionPPIs$DBORF[ConditionPPIs$ADORF==target])
}

