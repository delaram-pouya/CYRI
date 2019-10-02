## import needed functions and libraries
source('Code/Functions.R')
Initialize()

### Input data
ListOfConditions <- readRDS('Data/ListOfConditions_IS7.rds')
lapply(ListOfConditions, head)
GOterms <- read.delim('Data/GO/GO_terms_results.tsv', header = F)
colnames(GOterms) <- c('SGD_id', 'sys_name', 'Organism', 'stand_name', 'Description')


### you get very different results by counting or not counting the duplicates
combined_degrees <- lapply(ListOfConditions, .getComb_Degree)
lapply(combined_degrees, head)


combined_hub_thresholds <- unlist(lapply(combined_degrees, function(x) quantile(x$Freq, 0.95)))
ListOfHubs <- sapply(1:length(combined_degrees), function(i) 
  subset(combined_degrees[[i]], Freq>combined_hub_thresholds[i]), simplify = F)
names(ListOfHubs) <- names(ListOfConditions)
lapply(ListOfHubs, head)

combined_degree_Matrix <- getDegreeMatrix(combined_degrees)
compare_hubs <- getHubComparisonMatrix(ListOfHubs)

comb_result_matrix <- get_attributes(compare_hubs, combined_degree_Matrix) 

h2o2_increased <- subset(comb_result_matrix, Baseline==0 & H2O2 > 0 )
mms_increased <- subset(comb_result_matrix, Baseline==0 & MMS > 0 )
poorcarbon_increased <- subset(comb_result_matrix, Baseline==0 & PoorCarbon > 0 ) 

listOf_tfitTab <- readRDS('Data/listOf_tfitTab_DE_analysis.rds')
lapply(listOf_tfitTab, head)
listOf_DE <- lapply(listOf_tfitTab, function(x) subset(x, (logFC>2 | logFC< (-2)) & (P.Value<0.05) & (adj.P.Val<0.05)  ))
lapply(listOf_DE, head)


DE_genes_Base.vs.H2O2 <- getDEgenes(rownames(listOf_DE$Base.vs.H2O2)) #none
DE_genes_Base.vs.MMS <- getDEgenes(rownames(listOf_DE$Base.vs.MMS)) 
DE_genes_Base.vs.PoorCarbon <- getDEgenes(rownames(listOf_DE$Base.vs.PoorCarbon))


h2o2_increased[h2o2_increased$stand_name %in% DE_genes_Base.vs.H2O2, ] #none
mms_increased[mms_increased$stand_name %in% DE_genes_Base.vs.MMS, ] # none
poorcarbon_increased[poorcarbon_increased$stand_name %in% DE_genes_Base.vs.PoorCarbon, ] ## 2 are DE

### MMS
mms_increased$logFCDeg <- (mms_increased$Freq_MMS+0.1)/(mms_increased$Freq_Baseline+0.1)
mms_increased <- mms_increased[order(mms_increased$logFCDeg, decreasing = T),]

### H2O2
h2o2_increased$logFCDeg <- (h2o2_increased$Freq_H2O2+0.1)/(h2o2_increased$Freq_Baseline+0.1)
h2o2_increased <- h2o2_increased[order(h2o2_increased$logFCDeg, decreasing = T),]

### PoorCarbon
poorcarbon_increased$logFCDeg <- (poorcarbon_increased$Freq_PoorCarbon+0.1)/(poorcarbon_increased$Freq_Baseline+0.1)
poorcarbon_increased <- poorcarbon_increased[order(poorcarbon_increased$logFCDeg, decreasing = T),]



write.csv(mms_increased, 'Results/Result_Tables/Combined/mms_increased.csv')
write.csv(h2o2_increased, 'Results/Result_Tables/Combined/h2o2_increased.csv')
write.csv(poorcarbon_increased, 'Results/Result_Tables/Combined/poorcarbon_increased.csv')








