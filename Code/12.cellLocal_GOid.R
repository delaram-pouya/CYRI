## Possible Workflow
# 1- have a list of proteins which are present in both states, 
# 2- Find the set of interacting proteins with them in both conditions
# 3- Find the localization of each protein based on its interacting pairs in each condition 
# - Run enrichment analysis of these 2 sets → MMS-set , Baseline-set
# - Majority vote
# 4- for each protein, check if the cellular compartments are actually different in these 2 conditions

source('Code/Functions.R')
Initialize()

## Function
getListOfInteractingPairs <- function(sharedGenes , aConditionPPImat){
  sapply(1:nrow(sharedGenes), 
         function(i){
           geneOfInterest <- sharedGenes$ORF[i]
           interacting_pair_DBORF <- unique(aConditionPPImat$DBORF[aConditionPPImat$ADORF==geneOfInterest])
           interacting_pair_ADORF <- unique(aConditionPPImat$ADORF[aConditionPPImat$DBORF==geneOfInterest])
           interacting_pair <- union(interacting_pair_DBORF, interacting_pair_ADORF)
         }, simplify = F)
}


getListOfLocalizationForPairConditions <- function(stress_base_ORFS){
  sapply(1:nrow(stress_base_ORFS), function(i){
    
    genesToFindLocalization_baseline <- stress_base_ORFS$Interacting_pairs_baseline[[i]]
    genesToFindLocalization_stress <- stress_base_ORFS$Interacting_pairs_stress[[i]]
    
    baseline_localiz <- Gomap_localization$GOnumber[Gomap_localization$sys_name %in% genesToFindLocalization_baseline]
    stress_localiz <- Gomap_localization$GOnumber[Gomap_localization$sys_name %in% genesToFindLocalization_stress]

    baseline_localiz <- baseline_localiz[!is.na(baseline_localiz) & baseline_localiz!='' ]
    stress_localiz <- stress_localiz[!is.na(stress_localiz) & stress_localiz!='']
    
    list_localiz <- list(  round(table(baseline_localiz)/sum(table(baseline_localiz)), 2) , 
                         round(table(stress_localiz)/sum(table(stress_localiz)), 2)  )
    names(list_localiz) <- c('baseline_localiz', 'stress_localiz')
    
    return(list_localiz)
  }, simplify = F)
}

getHighestEnriched_baseline <- function(index, conditional_localization_pair){
  sample <- conditional_localization_pair[[index]]
  sample_baseline <- data.frame(sample$baseline_localiz)
  enrichedGoTerm <- ifelse(nrow(sample_baseline)>0 , 
                           as.character(sample_baseline[order(sample_baseline$Freq, decreasing = T), ][1,]$baseline_localiz),
                           '')
  enrichedGoTerm <- ifelse(length(enrichedGoTerm)>0 & !is.na(enrichedGoTerm) & enrichedGoTerm!='', enrichedGoTerm, '')
  return(enrichedGoTerm)
}

getHighestEnriched_stress <- function(index, conditional_localization_pair){
  sample <- conditional_localization_pair[[index]]
  sample_stress <- data.frame(sample$stress_localiz)
  enrichedGoTerm <- ifelse(nrow(sample_stress)>0 , 
                           as.character(sample_stress[order(sample_stress$Freq, decreasing = T), ][1,]$stress_localiz),
                           '')
  enrichedGoTerm <- ifelse(length(enrichedGoTerm)>0 & !is.na(enrichedGoTerm) & enrichedGoTerm!='', enrichedGoTerm, '')
  return(enrichedGoTerm)
}


getEnriched_GOterms <- function(index, conditional_localization_pair){
  
  ## enrichment -> freq higer in the new condition 0.3
  sample <- conditional_localization_pair[[index]]
  sample_baseline <- data.frame(sample$baseline_localiz)
  sample_stress <- data.frame(sample$stress_localiz)
  
  LOG_FOLD_CHANGE_CUTOFF = 2.3
  ## if the gene is also present in the baseline condition -> log2(stress/baseline >1)
  presentInBaseline <- as.character(sample_stress$stress_localiz[sample_stress$stress_localiz %in% 
                                                                   sample_baseline$baseline_localiz])
  foldChange <- log2((sample_stress$Freq[sample_stress$stress_localiz %in% presentInBaseline]+1e-4) /
                       (sample_baseline$Freq[sample_baseline$baseline_localiz %in% presentInBaseline]+1e-4))
  enriched_terms_present <- as.character(sample_stress$stress_localiz[sample_stress$stress_localiz %in% 
                                                                        presentInBaseline][foldChange>LOG_FOLD_CHANGE_CUTOFF])
  
  
  ## if not present -> freq > 0.3
  FREQUENCY_THRESHOLD = 0.7
  absentInBaseline <- as.character(sample_stress$stress_localiz[! sample_stress$stress_localiz %in% 
                                                                  sample_baseline$baseline_localiz])
  isAdequate <- absentInBaseline[sample_stress$Freq[sample_stress$stress_localiz %in% absentInBaseline]>FREQUENCY_THRESHOLD]
  enriched_terms_absent <- ifelse(length(isAdequate)>0, absentInBaseline,'')
  
  enrichedGO <- c(enriched_terms_absent, enriched_terms_present)
  enrichedGO <- enrichedGO[enrichedGO != '' &  !is.na(enrichedGO)]
  
  toReturn <- ifelse(length(enrichedGO)>0, enrichedGO, '')
  return(toReturn)
}


## get a list of proetins in each condition and their degrees
ListOfConditions <- readRDS('Data/ListOfConditions_IS7.rds')
combined_degrees <- lapply(ListOfConditions, .getComb_Degree)
lapply(combined_degrees, head)


###### GO slim mapping to map genes to pathways
GOmap <- read.delim('Data/GO/go_slim_mapping.tab', header = F)
colnames(GOmap)<- c('sys_name','stand_name','SGD_id', 'CFP','Pathway', 'GOnumber', 'type')
Gomap_localization <- GOmap[GOmap$CFP=='C', ]
table(Gomap_localization$GOnumber=='')

AD_degrees <- lapply(ListOfConditions, .getAD_Degrees)
DB_degrees <- lapply(ListOfConditions, .getDB_Degrees)
AD_Genes <- lapply(AD_degrees, function(x) data.frame(Gene=x$Gene))
DB_Genes <- lapply(DB_degrees, function(x) data.frame(Gene=x$Gene))
all_possible_genes <- data.frame(Genes=unique(c(as.character(unlist(AD_Genes)),
                                                as.character(unlist(DB_Genes)))))

### make this list for all genes
all_genes_localization <- merge(all_possible_genes, Gomap_localization, 
                                by.x='Genes', by.y='sys_name' , all.x=T, sort=F)

sum(is.na(all_genes_localization$GOnumber))
nrow(all_genes_localization)

# 1- have a list of proteins which are present in both states, 
mms_base_ORFS <- data.frame(ORF=intersect(combined_degrees$Baseline$Gene, combined_degrees$MMS$Gene))
h2o2_base_ORFS <- data.frame(ORF=intersect(combined_degrees$Baseline$Gene, combined_degrees$H2O2$Gene))
poorcarbon_base_ORFS <- data.frame(ORF=intersect(combined_degrees$Baseline$Gene, combined_degrees$PoorCarbon$Gene))

### MMS and baseline
mms_base_ORFS$Interacting_pairs_baseline <- getListOfInteractingPairs(mms_base_ORFS , ListOfConditions$Baseline)
mms_base_ORFS$Interacting_pairs_stress <- getListOfInteractingPairs(mms_base_ORFS , ListOfConditions$MMS)

### H2O2 and baseline
h2o2_base_ORFS$Interacting_pairs_baseline <- getListOfInteractingPairs(h2o2_base_ORFS , ListOfConditions$Baseline)
h2o2_base_ORFS$Interacting_pairs_stress <- getListOfInteractingPairs(h2o2_base_ORFS , ListOfConditions$H2O2)

### PoorCarbon and baseline
poorcarbon_base_ORFS$Interacting_pairs_baseline <- getListOfInteractingPairs(poorcarbon_base_ORFS , ListOfConditions$Baseline)
poorcarbon_base_ORFS$Interacting_pairs_stress <- getListOfInteractingPairs(poorcarbon_base_ORFS , ListOfConditions$PoorCarbon)

tail(mms_base_ORFS)
tail(h2o2_base_ORFS)
tail(poorcarbon_base_ORFS)

### how exactly should I compare these two list?
mms_base_ORFS_localizPair <- getListOfLocalizationForPairConditions(mms_base_ORFS)
h2o2_base_ORFS_localizPair <- getListOfLocalizationForPairConditions(h2o2_base_ORFS)
poorcarbon_base_localizPair <- getListOfLocalizationForPairConditions(poorcarbon_base_ORFS)


mms_base_ORFS$enriched <- sapply(1:length(mms_base_ORFS_localizPair),
                                 function(i){getEnriched_GOterms(i, mms_base_ORFS_localizPair)})
h2o2_base_ORFS$enriched <- sapply(1:length(h2o2_base_ORFS_localizPair),
                                  function(i){getEnriched_GOterms(i, h2o2_base_ORFS_localizPair)})
poorcarbon_base_ORFS$enriched <- sapply(1:length(poorcarbon_base_localizPair), 
                                        function(i){getEnriched_GOterms(i, poorcarbon_base_localizPair)})


mms_base_ORFS$numberOfBaseLine <- sapply(1:nrow(mms_base_ORFS),function(i) length(mms_base_ORFS$Interacting_pairs_baseline[[i]]))
mms_base_ORFS$numberOfStress <- sapply(1:nrow(mms_base_ORFS),function(i) length(mms_base_ORFS$Interacting_pairs_stress[[i]]))

h2o2_base_ORFS$numberOfBaseLine <- sapply(1:nrow(h2o2_base_ORFS),function(i) length(h2o2_base_ORFS$Interacting_pairs_baseline[[i]]))
h2o2_base_ORFS$numberOfStress <- sapply(1:nrow(h2o2_base_ORFS),function(i) length(h2o2_base_ORFS$Interacting_pairs_stress[[i]]))

poorcarbon_base_ORFS$numberOfBaseLine <- sapply(1:nrow(poorcarbon_base_ORFS),function(i) length(poorcarbon_base_ORFS$Interacting_pairs_baseline[[i]]))
poorcarbon_base_ORFS$numberOfStress <- sapply(1:nrow(poorcarbon_base_ORFS),function(i) length(poorcarbon_base_ORFS$Interacting_pairs_stress[[i]]))


# !colnames(mms_base_ORFS) %in%  c('Interacting_pairs_baseline', 'Interacting_pairs_stress')
mms_base_ORFS_2 <- mms_base_ORFS[ mms_base_ORFS$enriched != '' & mms_base_ORFS$numberOfBaseLine>2 & mms_base_ORFS$numberOfStress>2,]
h2o2_base_ORFS_2 <- h2o2_base_ORFS[h2o2_base_ORFS$enriched != '' & h2o2_base_ORFS$numberOfBaseLine>2 & h2o2_base_ORFS$numberOfStress>2,]
poorcarbon_base_ORFS_2 <- poorcarbon_base_ORFS[poorcarbon_base_ORFS$enriched != '' & poorcarbon_base_ORFS$numberOfBaseLine>2 & poorcarbon_base_ORFS$numberOfStress>2,]


### mapping GO terms to their pathways (C)

GO_pathway_dict <- unique(GOmap[,c('Pathway', 'GOnumber')])

mms_base_ORFS_2 <- merge(mms_base_ORFS_2, GO_pathway_dict, by.x=c('enriched'), by.y=c('GOnumber'), all.x=T)
h2o2_base_ORFS_2 <- merge(h2o2_base_ORFS_2, GO_pathway_dict, by.x=c('enriched'), by.y=c('GOnumber'), all.x=T)
poorcarbon_base_ORFS_2 <- merge(poorcarbon_base_ORFS_2, GO_pathway_dict, by.x=c('enriched'), by.y=c('GOnumber'), all.x=T)


nrow(mms_base_ORFS_2)
nrow(h2o2_base_ORFS_2)
nrow(poorcarbon_base_ORFS_2)
View(poorcarbon_base_ORFS_2)
View(h2o2_base_ORFS_2)
## -----top5 ---
## fisher exact test
## FuncAssiciate


### how to compre the 2 conditions??
### is the data enough?? -> how many are not getting mapped??
### what other enrichment softwares are available?



