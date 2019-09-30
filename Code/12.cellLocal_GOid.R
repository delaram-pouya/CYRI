## Possible Workflow
# 1- have a list of proteins which are present in both states, 
# 2- Find the set of interacting proteins with them in both conditions
# 3- Find the localization of each protein based on its interacting pairs in each condition 
# - Run enrichment analysis of these 2 sets → MMS-set , Baseline-set
# - Majority vote, fisher exact test
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
    
    list_localiz <- list(table(baseline_localiz) , table(stress_localiz))
    names(list_localiz) <- c('baseline_localiz', 'stress_localiz')
    
    return(list_localiz)
  }, simplify = F)
}


## read the links
# http://www.biostathandbook.com/fishers.html
# https://www.r-bloggers.com/contingency-tables-–-fisher’s-exact-test/
# https://www.r-bloggers.com/tag/fishers-exact-test/

getFisherPvalues <- function(df_merged){
  sapply(1:nrow(df_merged), 
         function(GO_index){
           GO_count_baseline <- df_merged$Freq.baseline[GO_index]
           otherGO_count_baseline <- sum(df_merged$Freq.baseline[-GO_index])
           GO_count_stress <- df_merged$Freq.stress[GO_index]
           otherGO_count_stress <- sum(df_merged$Freq.stress[-GO_index])
           
           fisher_table <- data.frame(matrix(data=c(GO_count_baseline,
                                                    GO_count_stress, 
                                                    otherGO_count_baseline, 
                                                    otherGO_count_stress), 
                                             nrow=2, ncol= 2,byrow = T), row.names = c('GOterm', 'Others'))
           
           colnames(fisher_table) <- c('baseline', 'stress')
           fisher_res <- fisher.test(fisher_table)
           return(fisher_res$p.value)
         })
}


## iterate over all the proteins 
## for each protein 
## for each GO term 
## construct the frequency table 
## run the fisher-exact test to check of that condition is enriched for that protein
## return the enriched condition for each protein
## return the frequency of enriched protein in the baseline and stress condition
## return the frequency of that GO term in the genome-wide analysis

getEnriched_GOterms <- function(stress_base_ORFS_localizPair){
  sapply(1:length(stress_base_ORFS_localizPair), 
         function(index){
           flag= FALSE
           #print(index)
           goList_base.df <- data.frame(stress_base_ORFS_localizPair[[index]]$baseline_localiz)
           goList_stress.df <- data.frame(stress_base_ORFS_localizPair[[index]]$stress_localiz)
           
           if(nrow(goList_stress.df)==0 & nrow(goList_base.df)!=0) {
             goList_merged = data.frame(GO_terms=goList_base.df$baseline_localiz,
                                        Freq.baseline= goList_base.df$Freq, 
                                        Freq.stress=rep(0, nrow(goList_base.df)))
             
           }else if(nrow(goList_base.df)==0 & nrow(goList_stress.df)!=0){
             goList_merged = data.frame(GO_terms=goList_stress.df$stress_localiz,
                                        Freq.baseline= rep(0, nrow(goList_stress.df)), 
                                        Freq.stress=goList_stress.df$Freq)
             
           }else if( nrow(goList_stress.df)==0 & nrow(goList_base.df)==0 ){
             flag = TRUE
           }else{
             goList_merged = merge(goList_base.df, goList_stress.df, by.x='baseline_localiz', by.y='stress_localiz', all.x=T, all.y=T)
           }
           
           if(flag) return('') else{
             goList_merged[is.na(goList_merged)] <- 0 
             colnames(goList_merged)<- c('GO_terms', 'Freq.baseline', 'Freq.stress')
             goList_merged$fisher_pval <- getFisherPvalues(goList_merged)
             enriched_terms <- ifelse(goList_merged$fisher_pval<0.1,  
                                      as.character(goList_merged$GO_terms[goList_merged$fisher_pval<0.1]), '')
             
             enriched_terms <- enriched_terms[enriched_terms!='']
             if(length(enriched_terms)>0 ) print(index)
             toReturn <- ifelse(length(enriched_terms)>0, enriched_terms, '')
             return(toReturn)
           }       
         })
}



getGOtermFreq_inBaseline <- function(stress_base_ORFS_localizPair, stress_base_ORFS, ORF_NAME, GO_NAME){
  
  target_orf_info = stress_base_ORFS_localizPair[[which(stress_base_ORFS$ORF==ORF_NAME)]]
  
  b = data.frame(target_orf_info$baseline_localiz)
  count_in_baseline = b$Freq[b$baseline_localiz==GO_NAME]
  if(length(count_in_baseline)==0) count_in_baseline = 0 
  freq_in_basline = count_in_baseline/sum(b$Freq)
  return(round(freq_in_basline,2))
}



getGOtermFreq_inStress <- function(stress_base_ORFS_localizPair, stress_base_ORFS, ORF_NAME, GO_NAME){
  
  target_orf_info = stress_base_ORFS_localizPair[[which(stress_base_ORFS$ORF==ORF_NAME)]]
  
  s = data.frame(target_orf_info$stress_localiz)
  count_in_stress = s$Freq[s$stress_localiz==GO_NAME]
  if(length(count_in_stress)==0) count_in_stress = 0 
  freq_in_stress = count_in_stress/sum(s$Freq)
  return(round(freq_in_stress, 2))
}




addFrequencyAttrib_GO <- function(Results, stress_base_ORFS_localizPair, stress_base_ORFS){
  ### mapping GO terms to their pathways (C)
  Results = merge(Results, GO_pathway_dict, by.x=c('enriched_term'), by.y=c('GOnumber'), all.x=T)
  
  ## add the frequency of that GO term in the baseline
  Results$GOfreqBaseline = sapply(1:nrow(Results), 
                                  function(i) getGOtermFreq_inBaseline(stress_base_ORFS_localizPair, 
                                                                       stress_base_ORFS, 
                                                                       Results$ORF[i], 
                                                                       Results$enriched_term[i] ))
  
  ## add the frequency of that GO term in the stress condition
  Results$GOfreqStress = sapply(1:nrow(Results), 
                                function(i) getGOtermFreq_inStress(stress_base_ORFS_localizPair, 
                                                                   stress_base_ORFS, 
                                                                   Results$ORF[i], 
                                                                   Results$enriched_term[i] ))
  
  ## get the global frequency of the GO term in all the genes that we have 
  Results <- merge(Results, GlobalGOFreq, by.x='enriched_term', by.y= 'GOterms', all.x=T)
  return(Results)
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

GO_pathway_dict <- unique(GOmap[,c('Pathway', 'GOnumber')])

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


GlobalGOFreq <-data.frame(round(table(all_genes_localization$GOnumber)/nrow(all_genes_localization), 2))
colnames(GlobalGOFreq) = c('GOterms', 'Global_Freq')



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


##### Peforming the enrichment analysis based in the GO terms
mms_base_ORFS$enriched_term <- getEnriched_GOterms(mms_base_ORFS_localizPair)
mmsRes <- subset(mms_base_ORFS[mms_base_ORFS$enriched_term != '',], select= c(ORF, enriched_term))

h2o2_base_ORFS$enriched_term <- getEnriched_GOterms(h2o2_base_ORFS_localizPair)
h2o2Res <-subset(h2o2_base_ORFS[h2o2_base_ORFS$enriched_term != '',] , select= c(ORF, enriched_term))

poorcarbon_base_ORFS$enriched_term <- getEnriched_GOterms(poorcarbon_base_localizPair)
poorcarbonRes <- subset(poorcarbon_base_ORFS[poorcarbon_base_ORFS$enriched_term != '',], select= c(ORF, enriched_term))


####### Adding all the frequency attributes
mmsRes = addFrequencyAttrib_GO(mmsRes, mms_base_ORFS_localizPair, mms_base_ORFS )
h2o2Res= addFrequencyAttrib_GO(h2o2Res, h2o2_base_ORFS_localizPair, h2o2_base_ORFS )
poorcarbonRes = addFrequencyAttrib_GO(poorcarbonRes, poorcarbon_base_localizPair, poorcarbon_base_ORFS )

head(mmsRes)








