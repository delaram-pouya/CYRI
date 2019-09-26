## import needed functions and libraries
source('Code/Functions.R')
Initialize()

# ListOfConditions <- list(Baseline, H2O2, MMS, PoorCarbon)
ListOfConditions <- readRDS('Data/ListOfConditions_IS7.rds')
lapply(ListOfConditions, head)

listOfPPIs <- lapply(ListOfConditions, function(x) x$PPI)
names(listOfPPIs) <- names(ListOfConditions)
all_PPI_uniq <- unique(unlist(listOfPPIs))

### Comparing the 'PRESENCE' of PPI in different networks
compare_PPIs <- data.frame(do.call(cbind, lapply(ListOfConditions, function(x) ifelse(all_PPI_uniq %in% x$PPI, 1, 0))))
rownames(compare_PPIs) <- all_PPI_uniq
compare_PPIs$sum <- rowSums(compare_PPIs)
tail(compare_PPIs)

p_ppi_compare <- ggplot(compare_PPIs, aes(sum))+
  geom_histogram(bins=7, color='black', fill='yellow', alpha= 0.5)+theme_bw()+xlab('#Common PPI pairs')
p_ppi_compare

### Venn diagram
lapply(listOfPPIs, head)
vennPlot_PPI <- venn.diagram(listOfPPIs , NULL, fill=rainbow(4),
                            alpha=c(0.4,0.4,0.4,0.4), 
                            cex =2.2, cat.fontface=1, 
                            category.names=names(test), main = 'PPI-pairs comparison', main.cex = 2)
grid.draw(vennPlot_PPI)

getPairPPIMatrix <- function(PPIs){
  PPIs_tmp <- sapply(1:length(PPIs), 
                        function(i) {
                          colnames(PPIs[[i]])[2] <- paste0('Freq','_', names(PPIs)[i]);PPIs[[i]]}, 
                        simplify = F)
  names(PPIs_tmp) <- names(PPIs)
  PPIs_Matrix <- Reduce(function(x, y) merge(x, y,by.x='PPI',by.y='PPI' ,all=TRUE), PPIs_tmp)
  PPIs_Matrix[is.na(PPIs_Matrix)] <- 0
  return(PPIs_Matrix)
}

listOfPPIs_df <- lapply(listOfPPIs, function(x) data.frame(table(x)))
listOfPPIs_df <- lapply(listOfPPIs_df, 
                        function(x) {
                          colnames(x)[1]<-'PPI'
                          x <- x[order(x$Freq, decreasing = T),]
                          return(x)} )

PPIs_df <- getPairPPIMatrix(listOfPPIs_df)
head(PPIs_df)
saveRDS( PPIs_df, 'Results/PPI_df.rds')

