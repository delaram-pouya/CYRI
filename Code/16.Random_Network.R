### it this script, we'll try to generate a random network and 
### then check athe accuracy of the localization prediction based on that

source('Code/Functions.R')
source('Code/Enrichment_Functions.R')
Initialize()

clean.tab <- function(x) {
  df = data.frame(table(x$GOnumber))
  df = df[-1,]
  colnames(df) = c('GOterm', 'Count')
  return(df)}



##### Import data
GOmap <- read.delim('Data/GO/go_slim_mapping.tab', header = F)
colnames(GOmap)<- c('sys_name','stand_name','SGD_id', 'CFP','Pathway', 'GOnumber', 'type')
Gomap_localization <- GOmap[GOmap$CFP=='C', ]
GO_pathway_dict <- unique(Gomap_localization[,c('Pathway', 'GOnumber')])

Proteins_in_each_condition.go.Count.table <- readRDS('Results/Proteins_in_each_condition.go.Count.table.rds')
lapply(Proteins_in_each_condition.go.Count.table, head)

## get a list of proetins in each condition and their degrees
ListOfConditions <- readRDS('Data/ListOfConditions_IS7.rds')
combined_degrees <- lapply(ListOfConditions, .getComb_Degree)

Baseline <- ListOfConditions$Baseline
listOfADs <- Baseline$ADORF
listOfDBs <- Baseline$DBORF
pred.score <- rep(NA, 100)


for(i in 1:100) {
  
  print(i)
  randNet <- data.frame(ADORF=sample(listOfADs), DBORF=sample(listOfDBs))
  randNet$PPI <- paste0(randNet$ADORF, '_', randNet$DBORF)
  rand_combined_degrees <- .getComb_Degree(randNet)
  
  ##### making the mapping between gene names and GO terms for the random network 
  AD_degrees.rand <-  .getAD_Degrees(randNet)
  DB_degrees.rand <-  .getDB_Degrees(randNet)
  AD_Genes.rand <- data.frame(Gene=AD_degrees.rand$Gene)
  DB_Genes.rand <- data.frame(Gene=DB_degrees.rand$Gene)
  
  Proteins_in_rand.df <- data.frame(Gene=unique(c(as.character(AD_Genes.rand$Gene), as.character(DB_Genes.rand$Gene))))
  Proteins_in_rand.go <- merge(Proteins_in_rand.df, Gomap_localization, by.x='Gene', by.y='sys_name' , all.x=T, sort=F)
  Proteins_in_rand.go.Count.table <- clean.tab( Proteins_in_rand.go)
  
  ##############  Predicting the localization for proteins in the random network 
  rand_ORFS <- getCompletePairingTable(rand_combined_degrees,  randNet)
  rand_ORFS <- getCompleteEnrichments(rand_ORFS, Proteins_in_rand.go.Count.table)
  rand_ORFS_localized <- getFinalSubset(rand_ORFS)
  pred.score[i] <- sum(rand_ORFS_localized$Consistent.loc >0)/nrow(rand_ORFS_localized)  # 0.06521739
  print(pred.score[i])
}


pred.score <- readRDS('Results/predScores.rds')
x = data.frame(score=c(pred.score, 0.2419355), category=c(rep('random',length(pred.score)), 'baseline'))
ggplot(x, aes(y=score, x=category))+geom_boxplot(aes(fill=category))+xlab('')




baseline1 <- subset(Baseline, select=c('ADORF', 'DBORF', 'Confidence'))
baseline1$score = baseline1$Confidence
colnames(baseline1) = c('Bait','Prey','HEKScore','AP-MS Score')
write.csv(baseline1, '../Desktop/cyto_test.csv', quote = F, row.names = F)


### check this for degree controlled random network generation (is the current method label swapping?)
# https://dnac.ssri.duke.edu/r-labs/2016/08_random_graphs.php
# install cytoscape -> visualize -> color each node -> color with respect to visualization -> use official name
## check the ORFs which have been labeled as consisitent in the base but are not consistent in the random network
#### what was their neighboring proteins which led to this accurate prediction 

## plot the degree distribution for both the baseline and the random network -> both with/without duplicates

## for localization change prediction -> add the other pvalues as well 

  
  