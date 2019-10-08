
source('Code/Functions.R')
source('Code/Enrichment_Functions.R')
Initialize()


##### Import data
GOmap <- read.delim('Data/GO/go_slim_mapping.tab', header = F)
colnames(GOmap)<- c('sys_name','stand_name','SGD_id', 'CFP','Pathway', 'GOnumber', 'type')
Gomap_localization <- GOmap[GOmap$CFP=='C', ]
GO_pathway_dict <- unique(Gomap_localization[,c('Pathway', 'GOnumber')])

Proteins_in_each_condition.go.Count.table <- readRDS('Results/Proteins_in_each_condition.go.Count.table.rds')
lapply(Proteins_in_each_condition.go.Count.table, head)

## get a list of proetins in each condition and their degrees
ListOfConditions <- readRDS('Data/ListOfConditions_IS7.rds')
lapply(ListOfConditions, head)
combined_degrees <- lapply(ListOfConditions, .getComb_Degree)
lapply(combined_degrees, head)



############## Baseline
## make a list of ORFs in a condition & find their interacting pairs, map the protein names to the GO number
base_ORFS <- getCompletePairingTable(combined_degrees$Baseline,  ListOfConditions$Baseline)
### for all the ORFs, find the enriched GO term and its p.value 
base_ORFS <- getCompleteEnrichments(base_ORFS, Proteins_in_each_condition.go.Count.table$Baseline)
## subset of the protein of which a GO term has been enriched & find the true localization for the ORF accouring to GO dataset
## check if the enriched GO term based on the interacting pairs is consistent with the true localization
base_ORFS_localized <- getFinalSubset(base_ORFS)
head(base_ORFS_localized)
sum(base_ORFS_localized$Consistent.loc >0)/nrow(base_ORFS_localized)  # 0.2419355


############## MMS
mms_ORFS <- getCompletePairingTable(combined_degrees$MMS,  ListOfConditions$MMS)
mms_ORFS <- getCompleteEnrichments(mms_ORFS, Proteins_in_each_condition.go.Count.table$MMS)
mms_ORFS_localized <- getFinalSubset(mms_ORFS)
sum(mms_ORFS_localized$Consistent.loc >0)/nrow(mms_ORFS_localized) # 0.2285115


############## H2O2
h2o2_ORFS <- getCompletePairingTable(combined_degrees$H2O2,  ListOfConditions$H2O2)
h2o2_ORFS <- getCompleteEnrichments(h2o2_ORFS, Proteins_in_each_condition.go.Count.table$H2O2)
h2o2_ORFS_localized <- getFinalSubset(h2o2_ORFS)
sum(h2o2_ORFS_localized$Consistent.loc >0)/nrow(h2o2_ORFS_localized) # 0.2121896


############## PoorCarbon
poorcarbon_ORFS <- getCompletePairingTable(combined_degrees$PoorCarbon,  ListOfConditions$PoorCarbon)
poorcarbon_ORFS <- getCompleteEnrichments(poorcarbon_ORFS, Proteins_in_each_condition.go.Count.table$PoorCarbon)
poorcarbon_ORFS_localized <- getFinalSubset(poorcarbon_ORFS)
sum(poorcarbon_ORFS_localized$Consistent.loc >0)/nrow(poorcarbon_ORFS_localized) # 0.1921708



listOfGlobal_GO_enrichments = list(base_ORFS_localized, mms_ORFS_localized, h2o2_ORFS_localized, poorcarbon_ORFS_localized)
names(listOfGlobal_GO_enrichments) <- c('baseline', 'mms', 'h2o2', 'poorcarbon')
saveRDS(listOfGlobal_GO_enrichments, 'Results/list_Of_GO_enriched_Global_localizations.rds')



# calculate 2 other p-values for comparing the GO (or the protein??) with the global amount 
## add these 2 p.values to the table made in the other script for mms at least

## use this randomized network -> not label swapping 
# degree controlled randomization -> 
## odd ration between the protein 



####### visualization 

library(igraph)
set.seed(1)
gs <- list()
for (x in seq_len(20L)) {
  gs[[x]] <- erdos.renyi.game(sample(1:100, 1), p.or.m = runif(1))
  E(gs[[x]])$weight <- sample(1:5, ecount(gs[[x]]), T)
}
plot(gs[[1]], edge.width = E(gs[[1]])$weight) # plot first graph

