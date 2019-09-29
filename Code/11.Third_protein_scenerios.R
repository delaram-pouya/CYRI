
source('Code/Functions.R')
Initialize()

###### WORKFLOW
# check the PPIs which are present in both conditions -> 
# find the ones which in the condition, the confidence has increased
## find any possible third gene that has interaction with both the proteins
## check if this third gene's expression has increased in this new condition

## go through literature to find any biological relevance

## Functions:
## PTM 

### Is there any third protein that is interacting with poth pairs?
get3rdProteins <- function(PPI.highConf, targetCondition){
  ## for  (i in rows):
  sapply(1:nrow(PPI.highConf), 
         function(i){
           pair <- str_split_fixed(PPI.highConf$PPI_Baseline, '_', 2)[i, ]
           ## find the list of protein which are interacting with the first one
           interact_1 <- c(targetCondition$DBORF[targetCondition$ADORF==pair[1]], 
                           targetCondition$ADORF[targetCondition$DBORF==pair[1]])
           
           ## find the list of protein which are interacting with the second one
           interact_2 <- c(targetCondition$DBORF[targetCondition$ADORF==pair[2]], 
                           targetCondition$ADORF[targetCondition$DBORF==pair[2]])
           
           ## check if the intersecting proteins 
           # -> if intersecting with either of the pairs is ok -> just take  a union
           interact_both <- intersect(interact_1, interact_2)
           ToReturn <- ifelse(length(interact_both)>0 & pair[1]!=pair[2], interact_both, '' )
           return(ToReturn)
         }, simplify = T)
}


getCleanedTable <- function(PPI_highConf){
  ThirdProtein <- data.frame(table(PPI_highConf$ThirdProtein))
  ThirdProtein <- ThirdProtein[-1, ]
  return(ThirdProtein)  
}


expression <- read.csv('Data/genes.fpkm_tracking_clean.csv')
all_possible_genes_filled <- readRDS('mapGeneExp/all_possible_genes_exp.rds')
unmaped_indices <- readRDS('mapGeneExp/listOfUnmapped.rds')
ppi_exp_merged <- readRDS('Data/ppi_exp_merged.rds')

### Input data
ListOfConditions <- readRDS('Data/ListOfConditions_IS7.rds')
ListOfConditions2 <- lapply(ListOfConditions, function(x){subset(x, select= c('PPI', 'Confidence'))})
names(ListOfConditions2) <- names(ListOfConditions)

ListOfConditions2 <- sapply(1:length(ListOfConditions2), 
       function(i) {
         colnames(ListOfConditions2[[i]])=paste0(colnames(ListOfConditions2[[i]]),'_', names(ListOfConditions2)[i] )
         return(ListOfConditions2[[i]])}, 
       simplify = F)

names(ListOfConditions2) <- names(ListOfConditions)
lapply(ListOfConditions2, head)

ListOfConditions2$Baseline$Confidence_Baseline <- as.numeric(ListOfConditions2$Baseline$Confidence_Baseline)


###### DO WE NEED TO CONSIDER THE ONES WHICH ARE NOT PRESENT IN BOTH AS WELL??? -> zero instead of NA ??

PPI.base.mms <- merge(ListOfConditions2$Baseline, ListOfConditions2$MMS, by.x='PPI_Baseline', by.y='PPI_MMS')
PPI.base.poorcarbon <- merge(ListOfConditions2$Baseline, ListOfConditions2$PoorCarbon, by.x='PPI_Baseline', by.y='PPI_PoorCarbon')
PPI.base.h2o2 <- merge(ListOfConditions2$Baseline, ListOfConditions2$H2O2, by.x='PPI_Baseline', by.y='PPI_H2O2')

PPI.base.mms$diff.log2 <- log2(PPI.base.mms$Confidence_MMS/PPI.base.mms$Confidence_Baseline)
PPI.base.h2o2$diff.log2 <- log2(PPI.base.h2o2$Confidence_H2O2/PPI.base.h2o2$Confidence_Baseline)
PPI.base.poorcarbon$diff.log2 <- log2(PPI.base.poorcarbon$Confidence_PoorCarbon/PPI.base.poorcarbon$Confidence_Baseline)

hist(PPI.base.mms$diff.log2)
hist(PPI.base.poorcarbon$diff.log2)
hist(PPI.base.h2o2$diff.log2)


CONFIDENCE_CUT_OFF = 0.85

## list of proteins which in the stress condition, the confidence has increased
PPI.base.mms_highConf <- PPI.base.mms[PPI.base.mms$diff.log2 > quantile(PPI.base.mms$diff.log2, CONFIDENCE_CUT_OFF) & 
                                        PPI.base.mms$diff.log2 > 0,]
PPI.base.poorcarbon_highConf <- PPI.base.poorcarbon[PPI.base.poorcarbon$diff.log2 > quantile(PPI.base.poorcarbon$diff.log2, CONFIDENCE_CUT_OFF) & 
                                                      PPI.base.poorcarbon$diff.log2 > 0,]
PPI.base.h2o2_highConf <- PPI.base.h2o2[PPI.base.h2o2$diff.log2 > quantile(PPI.base.h2o2$diff.log2, CONFIDENCE_CUT_OFF) & 
                                          PPI.base.h2o2$diff.log2 > 0,]

PPI.base.mms_highConf$ThirdProtein <- get3rdProteins(PPI.base.mms_highConf, ListOfConditions$MMS)
PPI.base.poorcarbon_highConf$ThirdProtein <- get3rdProteins(PPI.base.poorcarbon_highConf, ListOfConditions$PoorCarbon)
PPI.base.h2o2_highConf$ThirdProtein <- get3rdProteins(PPI.base.h2o2_highConf, ListOfConditions$H2O2)



## for these candidates, check for the cases in which the 
##  Third protein is not '' & isn't the same as one of the proteins

isThirdProteinValid <- function(PPI.base.stress_highConf){
  cond1 <- ! PPI.base.stress_highConf$ThirdProtein %in% str_split_fixed(PPI.base.stress_highConf$PPI_Baseline, '_', 2)
  cond2 <- PPI.base.stress_highConf$ThirdProtein != ''
  return(cond1 & cond2)
}

## this is a list of third proteins which have interaction with both the protein 
# pair stress condition whose interaction has increased in the stress condition

PPI.base.mms_highConf <- PPI.base.mms_highConf[isThirdProteinValid(PPI.base.mms_highConf), ]
PPI.base.h2o2_highConf <- PPI.base.h2o2_highConf[isThirdProteinValid(PPI.base.h2o2_highConf), ]
PPI.base.poorcarbon_highConf <- PPI.base.poorcarbon_highConf[isThirdProteinValid(PPI.base.poorcarbon_highConf), ]


##### now we have to check if the expression level of this 3rd protein has actually changed
######  and then look for the increasing interactions

EXPRESSION_CUT_OFF = 0.2

### MMS 
PPI.base.mms_highConf2 <- merge(PPI.base.mms_highConf, ppi_exp_merged, by.x='ThirdProtein', by.y='systematic_name', all.x=T)
PPI.base.mms_highConf2 <- subset(PPI.base.mms_highConf2, 
                                 select=c(ThirdProtein,standard_name,  PPI_Baseline , Confidence_Baseline,  
                                          Confidence_MMS,diff.log2, exp_Baseline, exp_MMS,Exp_mms_base))
ThirdPr_final_mms <- subset(PPI.base.mms_highConf2, Exp_mms_base > EXPRESSION_CUT_OFF)
colnames(ThirdPr_final_mms) <- c('Third_protein_sys_name', 'standard_name', 'PPI_pair','Confidence_score_Baseline' , 
                                 'Confidence_score_MMS','conf.change.log2' , 'expression_Baseline'  , 'expression_MMS',  'exp.change.log2')
write.csv(ThirdPr_final_mms, '../Desktop/ThirdPr_final_mms.csv')

#### H2O2 --> Nothing :(
PPI.base.h2o2_highConf2 <- merge(PPI.base.h2o2_highConf, ppi_exp_merged, by.x='ThirdProtein', by.y='systematic_name', all.x=T)
PPI.base.h2o2_highConf2 <- subset(PPI.base.h2o2_highConf2, 
                                  select=c(ThirdProtein,standard_name,  PPI_Baseline , Confidence_Baseline, 
                                           Confidence_H2O2, diff.log2, exp_Baseline, exp_H2O2, Exp_h2o2_base))
ThirdPr_final_h2o2 <- subset(PPI.base.h2o2_highConf2, Exp_h2o2_base > EXPRESSION_CUT_OFF)
colnames(ThirdPr_final_h2o2) <- c('Third_protein_sys_name', 'standard_name', 'PPI_pair','Confidence_score_Baseline' , 
                                 'Confidence_score_H2O2','conf.change.log2' , 'expression_Baseline'  , 'expression_H2O2',  'exp.change.log2')



### Poorcarbon 
PPI.base.poorcarbon_highConf2 <- merge(PPI.base.poorcarbon_highConf, ppi_exp_merged, by.x='ThirdProtein', by.y='systematic_name', all.x=T)
PPI.base.poorcarbon_highConf2 <- subset(PPI.base.poorcarbon_highConf2, 
                                        select=c(ThirdProtein,standard_name,  PPI_Baseline , Confidence_Baseline,Confidence_PoorCarbon,
                                                 diff.log2, exp_Baseline, exp_PoorCarbon,Exp_poorcarbon_base))
ThirdPr_final_poorcarbon <-subset(PPI.base.poorcarbon_highConf2, Exp_poorcarbon_base > EXPRESSION_CUT_OFF)
colnames(ThirdPr_final_poorcarbon) <- c('Third_protein_sys_name', 'standard_name', 'PPI_pair','Confidence_score_Baseline' , 
                                  'Confidence_score_PoorCarbon','conf.change.log2' , 'expression_Baseline'  , 'expression_PoorCarbon',  'exp.change.log2')
write.csv(ThirdPr_final_poorcarbon, '../Desktop/ThirdPr_final_poorcarbon.csv')



#####################
ppi_exp_merged <- readRDS('ppi_exp_merged.rds')
ThirdProtein.mms <- getCleanedTable(PPI.base.mms_highConf)
ThirdProtein.mms <- merge(ThirdProtein.mms, ppi_exp_merged, by.x='Var1', by.y='systematic_name', all.x=T)

## check for each condition, which protein's expression has increased 
##    -> list of increased expression in mms/baseline expression 

## expression of first 2 genes don't chnage
### p values

ThirdPr_final_mms
ThirdPr_final_h2o2
ThirdPr_final_poorcarbon

### checking the Pvalues for the candidates
expression[expression$gene_short_name == ThirdPr_final_mms$standard_name,]$MMS_P_value
expression[expression$gene_short_name %in% ThirdPr_final_poorcarbon$Third_protein_sys_name,]
all_possible_genes <- readRDS('Data/all_possible_genes_exp.rds')  ## coudn't find the Pvalue for the poor carbon 

indices <- unlist(lapply(str_split(expression$gene_short_name, ',' ), function(x){ sum(ThirdPr_final_poorcarbon$standard_name%in% x)} ))
expression[indices, ]$Poor.Carbon_P_value


