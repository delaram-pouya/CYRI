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


index = 200
length(mms_base_ORFS_localizPair)
nrow(mms_base_ORFS)
mms_base_ORFS_localizPair[[index]]

sapply(1:length(mms_base_ORFS_localizPair), 
       function(index){
         flag= FALSE
         print(index)
         goList_base.df <- data.frame(mms_base_ORFS_localizPair[[index]]$baseline_localiz)
         goList_stress.df <- data.frame(mms_base_ORFS_localizPair[[index]]$stress_localiz)
         
         # if(nrow(goList_stress.df)==0) goList_stress.df=getNAdf() 
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
           enriched_terms <- ifelse(goList_merged$fisher_pval<0.05,  
                                    goList_merged$GO_terms[goList_merged$fisher_pval<0.05], '')
           
           enriched_terms <- enriched_terms[enriched_terms!='']
           toReturn <- ifelse(length(enriched_terms)>0, enriched_terms, '')
           return(toReturn)
         }       
       })


getNAdf <- function(){ data.frame(stress_localiz=NA, Freq=NA)}

getEnriched_GOterms <- function(){
}




### Do the fisher test here
## iterate over all the proteins 
## for each protein 
  ## for each GO term 
    ## construct the frequency table 
    ## run the fisher-exact test to check of that condition is enriched for that protein
    ## return the enriched condition for each protein
    ## return the frequency of enriched protein in the baseline and stress condition
    ## return the frequency of that GO term in the genome-wide analysis




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



