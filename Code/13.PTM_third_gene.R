# https://yeastkinome.org/supplementary.php

source('Code/Functions.R')
Initialize()

## Possible workflow:
# 1. Find all the possible third genes which expression is changing + is a PTM protein
# 2. Check all the hub proteins (not only hubs in this case) that those candidates interact with
# 3. Check the literature if you can find anything meaningfull 
# 4. Check the PPI-pairs for that hub -> are there any new interactions formed in the new condition?


## function
isInDifExpList <- function(difExp_geneNames){
  ifelse(  (all_possible_genes_filled$Genes %in% unlist(str_split(difExp_geneNames, ',')) & 
              !is.na(all_possible_genes_filled$Genes) &  all_possible_genes_filled$Genes != '') |
             
             (all_possible_genes_filled$stand_name %in% unlist(str_split(difExp_geneNames, ',')) &
                !is.na(all_possible_genes_filled$stand_name) &  all_possible_genes_filled$stand_name != '') , 
           
           TRUE, FALSE)
}

ListOfConditions <- readRDS('Data/ListOfConditions_IS7.rds')
combined_degrees <- lapply(ListOfConditions, .getComb_Degree)
lapply(ListOfConditions, head)
lapply(combined_degrees, head)


kinaseDB <- read.csv('Data/kinase_phospho_table.csv')
AllKinaseName <- c(kinaseDB$Orf.Name, kinaseDB$Gene.Name)
table(kinaseDB$Category)

### find the gene expressions which are changing -> possible 3rd genes
expression<- read.csv('Data/genes.fpkm_tracking_clean.csv')
all_possible_genes_filled <- readRDS('mapGeneExp/all_possible_genes_exp.rds')
all_possible_genes_filled[all_possible_genes_filled==-1] <- NA
all_possible_genes_filled[all_possible_genes_filled=='#DIV/0!'] <- 1
head(all_possible_genes_filled)



##### Thresholds for DE analysis
FOLD_CHANGE_THRESHOLD = 1.5
P_VAL_THRESHOLD = 0.2

listOf_tfitTab <- readRDS('Data/listOf_tfitTab_DE_analysis.rds')
listOf_tfitTab <- listOf_tfitTab[1:3]
listOf_tfitTab_difExp <- lapply(listOf_tfitTab, 
                                function(x) subset(x, ( logFC< -FOLD_CHANGE_THRESHOLD | logFC > FOLD_CHANGE_THRESHOLD) & 
                                                     P.Value<P_VAL_THRESHOLD))

difExp_trackIds <- lapply(listOf_tfitTab_difExp, rownames)
names(difExp_trackIds)
difExp_geneNames <- sapply(1:3, function(i)expression$gene_short_name[expression$tracking_id %in% difExp_trackIds[[i]]])
names(difExp_geneNames) <- names(difExp_trackIds)
lapply(difExp_geneNames, length)


all_possible_genes_filled$MMS_difExp <- isInDifExpList(difExp_geneNames$Base.vs.MMS)
all_possible_genes_filled$H2O2_difExp <- isInDifExpList(difExp_geneNames$Base.vs.H2O2)
all_possible_genes_filled$PoorCarbon_difExp <- isInDifExpList(difExp_geneNames$Base.vs.PoorCarbon)

MMS_difExp_df <- subset(all_possible_genes_filled, MMS_difExp)
H2O2_difExp_df <- subset(all_possible_genes_filled, H2O2_difExp)
PoorCarbon_difExp_df <- subset(all_possible_genes_filled, PoorCarbon_difExp)

#  all the possible third genes which expression is changing AND is a PTM protein
PTM_genes_MMS <- MMS_difExp_df$Genes[MMS_difExp_df$Genes %in% AllKinaseName]  # "YGL158W"
PTM_genes_H2O2 <- H2O2_difExp_df$Genes[H2O2_difExp_df$Genes %in% AllKinaseName] # nothing
PTM_genes_PoorCarbon <- PoorCarbon_difExp_df$Genes[PoorCarbon_difExp_df$Genes %in% AllKinaseName]



#### Finding pairs of those PTM proteins found in the previous section
############ MMS --> "YGL158W"
PTM_MMS_baseline_pairs <- c(ListOfConditions$MMS$DBORF[ListOfConditions$MMS$ADORF==PTM_genes_MMS], 
                            ListOfConditions$MMS$ADORF[ListOfConditions$MMS$DBORF==PTM_genes_MMS])

PTM_MMS_Stress_pair <- c(ListOfConditions$Baseline$DBORF[ListOfConditions$Baseline$ADORF==PTM_genes_MMS], 
                         ListOfConditions$Baseline$ADORF[ListOfConditions$Baseline$DBORF==PTM_genes_MMS])



## The assumption is that the PTM protein's expression has changed in the Stress condition and
### by this change -> new proteins are getting modified in the stress condition which weren't getting 
### modified before in the baseline condition -> we need to check if this newly affected protein interaction has 
### by getting modified in the Stress condition

ptm_effected_protein_mms <- PTM_MMS_Stress_pair[!PTM_MMS_Stress_pair %in% PTM_MMS_baseline_pairs]

ptm_df <- data.frame(PTM=rep(PTM_genes_MMS,2), ptm_effected=ptm_effected_protein_mms)
### Final answare for the PTM affected proteins in the new condition
b = combined_degrees$Baseline[combined_degrees$Baseline$Gene %in% ptm_effected_protein_mms, ]
colnames(b)[2] = 'deg_base'
pc = combined_degrees$MMS[combined_degrees$MMS$Gene %in% ptm_effected_protein_mms, ]
colnames(pc)[2] = 'deg_stress'
ptm_df_1 <- merge(ptm_df, b, by.x='ptm_effected', by.y='Gene', all.x=T)
ptm_df_mms <- merge(ptm_df_1, pc, by.x='ptm_effected', by.y='Gene', all.x=T)
ptm_df_mms[is.na(ptm_df_mms)] <- 0
View(ptm_df_mms)



#### Poorcarbon

ptm_effected_protein_poorcarbon <- sapply(1:length(PTM_genes_PoorCarbon), 
       function(i){
         PTM_gene <- PTM_genes_PoorCarbon[i]
         stress_pairs <- c(ListOfConditions$PoorCarbon$DBORF[ListOfConditions$PoorCarbon$ADORF==PTM_gene], 
           ListOfConditions$PoorCarbon$ADORF[ListOfConditions$PoorCarbon$DBORF==PTM_gene])
         
         baseline_pairs <- c(ListOfConditions$Baseline$DBORF[ListOfConditions$Baseline$ADORF==PTM_gene], 
           ListOfConditions$Baseline$ADORF[ListOfConditions$Baseline$DBORF==PTM_gene])
         
         
         ## we want the ones which the PTM-protein interacts with 
         ##   in the stress condition and not in the baseline
         stress_effected_protein <- stress_pairs[!stress_pairs %in% baseline_pairs]
         ifelse(length(stress_effected_protein)>0, stress_effected_protein, '' )
         }, simplify = T)


ptm_poorcarbon <- PTM_genes_PoorCarbon[ptm_effected_protein_poorcarbon != '']  ## PTM proteins
ptm_effected_protein_poorcarbon <- ptm_effected_protein_poorcarbon[ptm_effected_protein_poorcarbon != '']

ptm_df <- data.frame(PTM=ptm_poorcarbon, ptm_effected=ptm_effected_protein_poorcarbon) ## check literature


### Final answare for the PTM affected proteins in the new condition
b = combined_degrees$Baseline[combined_degrees$Baseline$Gene %in% ptm_effected_protein_poorcarbon, ]
colnames(b)[2] = 'deg_base'
pc = combined_degrees$PoorCarbon[combined_degrees$PoorCarbon$Gene %in% ptm_effected_protein_poorcarbon, ]
colnames(pc)[2] = 'deg_stress'
ptm_df_1 <- merge(ptm_df, b, by.x='ptm_effected', by.y='Gene', all.x=T)
ptm_df_poorcarbon <- merge(ptm_df_1, pc, by.x='ptm_effected', by.y='Gene', all.x=T)
ptm_df_poorcarbon[is.na(ptm_df_poorcarbon)] <- 0
View(ptm_df_poorcarbon)

ptm_df_poorcarbon$PTM
head(kinaseDB)
merge(ptm_df_poorcarbon, kinaseDB, by.x='PTM', by.y= 'Orf.Name')
merge(ptm_df_mms, kinaseDB, by.x='PTM', by.y= 'Orf.Name')


#### Poorcarbon
getInteractingPairs(ptm_df_poorcarbon$ptm_effected[3],ListOfConditions$Baseline )
getInteractingPairs(ptm_df_poorcarbon$ptm_effected[3],ListOfConditions$PoorCarbon )

### MMS
getInteractingPairs(ptm_df_mms$ptm_effected[2],ListOfConditions$Baseline )
getInteractingPairs(ptm_df_mms$ptm_effected[2],ListOfConditions$MMS )

## Tonight:
# read and document this tonight -> write out all the pairs -> check the literature





