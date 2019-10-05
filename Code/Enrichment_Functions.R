
getListOfInteractingPairs <- function(sharedGenes , aConditionPPImat){
  sapply(1:nrow(sharedGenes), 
         function(i){
           geneOfInterest <- sharedGenes$ORF[i]
           interacting_pair_DBORF <- unique(aConditionPPImat$DBORF[aConditionPPImat$ADORF==geneOfInterest])
           interacting_pair_ADORF <- unique(aConditionPPImat$ADORF[aConditionPPImat$DBORF==geneOfInterest])
           interacting_pair <- union(interacting_pair_DBORF, interacting_pair_ADORF)
         }, simplify = F)
}


getListOfLocalizationForSingleCondition <- function(ORF_list){
  sapply(1:nrow(ORF_list), function(i){
    genesToFindLocalization <- ORF_list$Interacting_pairs[[i]]
    localiz <- Gomap_localization$GOnumber[Gomap_localization$sys_name %in% genesToFindLocalization]
    localiz <- localiz[!is.na(localiz) & localiz!='' ]
    return(localiz)
  }, simplify = F)
}


### for all the GOterms associated with that ORF
getFisherPvalues_globalLocalization <- function(proteinNiche_list, Global.go.count.table ){
  sapply(1:nrow(proteinNiche_list), function(i){
    GO_term <- as.character(proteinNiche_list$GO_term[i])
    
    protein_niche_index = which(proteinNiche_list$GO_term==GO_term)
    GO_count_proteinNiche <- proteinNiche_list$Freq[ protein_niche_index ]
    GO_count_other <- sum(proteinNiche_list$Freq[-protein_niche_index])
    
    global_go_index = which(Global.go.count.table$GOterm == GO_term)
    global_GO_count <- Global.go.count.table$Count[ global_go_index]
    global_other_count <- sum(Global.go.count.table$Count[ -global_go_index])
    
    fisher_table = data.frame(matrix(c(GO_count_proteinNiche,
                                       global_GO_count,
                                       GO_count_other,
                                       global_other_count),
                                     nrow=2, ncol=2,byrow = T), row.names = c('GOterm', 'Others'))
    colnames(fisher_table) <- c('protein_niche', 'Global')
    fisher_res = fisher.test(fisher_table)
    return(fisher_res$p.value)
  }, simplify = T)
}



getEnriched_GOterms.global <- function(i, ORFS, Global.go.count.table){
  
  proteinNiche_list <- data.frame(table(ORFS$Interacting_pairs.localiz[[i]]))
  colnames(proteinNiche_list)[1] = c('GO_term')
  
  P_VALUE_THRESHOLD = 0.05
  pValuesForGOs <- as.numeric(getFisherPvalues_globalLocalization(proteinNiche_list, 
                                                                  Global.go.count.table))
  ### Find the enriched GO term for a protein
  enriched_GOs <- as.character(proteinNiche_list$GO_term[pValuesForGOs < P_VALUE_THRESHOLD])
  Pvalues <- pValuesForGOs[pValuesForGOs < P_VALUE_THRESHOLD]
  toReturn <- list(enriched_GOs, Pvalues)
  names(toReturn) <- c('GO', 'p.value')
  return(toReturn)}




getTrueLocalizations <- function(localized_ORFs){
  sapply(1:nrow(localized_ORFs), 
         function(i){
           an.ORF = localized_ORFs$ORF[i]
           trueLocalization = Gomap_localization$GOnumber[Gomap_localization$sys_name == an.ORF]
           return(trueLocalization)}, simplify = F)
}


isLocalizationCorrect <- function(index, ORFS_localized){
  ORF.info <-ORFS_localized[index, ]
  sum(unlist(ORF.info$enriched_GO) %in% unlist(ORF.info$true.localization))
}


getCompletePairingTable <- function(combined_degree.tab,aCondition.tab){
  
  ## make a list of ORFs in a condition
  ## find their interacting pairs, map the protein names to the GO number
  ORFS <- data.frame(ORF=combined_degree.tab$Gene)
  ORFS$Interacting_pairs <- getListOfInteractingPairs(ORFS , aCondition.tab)
  ORFS$Interacting_pairs.localiz <- getListOfLocalizationForSingleCondition(ORFS)
  return(ORFS)
  
}


getCompleteEnrichments <- function(ORFS, a.condition.GO.table ){
  ### for all the ORFs, find the enriched GO term and its p.value 
  enriched_pVal <- sapply(1:nrow(ORFS), function(i)  getEnriched_GOterms.global(i, ORFS, a.condition.GO.table), simplify = F)
  
  ORFS$enriched_GO <- lapply(enriched_pVal, function(x) x$GO)
  ORFS$p.values.global <- lapply(enriched_pVal, function(x) x$p.value)
  isEmpty <- unlist(lapply(enriched_pVal, function(x) ifelse(length(x$GO)==0, T, F)))
  ORFS$enriched_GO[isEmpty] <- ''
  ORFS$p.values.global[isEmpty]  <- ''
  return(ORFS)
}



getFinalSubset <- function(ORFS){
  ## subset of the protein of which a GO term has been enriched
  ORFS_localized <- ORFS[ORFS$enriched_GO != '', c('ORF', 'enriched_GO' , 'p.values.global')]
  
  ## find the true localization for the ORF accouring to GO dataset
  ## check if the enriched GO term based on the interacting pairs is 
  ## consistent with the true localization
  ORFS_localized$true.localization <- getTrueLocalizations(ORFS_localized)
  ORFS_localized$Consistent.loc <- sapply(1:nrow(ORFS_localized), function(i) isLocalizationCorrect(i, ORFS_localized))
  return(ORFS_localized)
}

