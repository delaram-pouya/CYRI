
## in this script, we'll try to compare the localization data provided in 
## https://www.nature.com/articles/ncb2549#s1 paper by brown lab, 
## with the localization changes by the MMS compound for the DNA damage stress condition
## the supplemetary table 2 is being used here

source('Code/Functions.R')
Initialize()

### Functions
get.Stress.Base.freq.table <- function(i){
  print(MMS.DataBrown_included[MMS.DataBrown_included$Systematic.ORF == checkWhyNotEnriched[i,]$ORF,])
  cat("\n")
  
  genesToMapToGO.base <- unlist(checkWhyNotEnriched[i, ]$Interacting_pairs_baseline)
  df_base = data.frame(table(Gomap_localization$GOnumber[Gomap_localization$sys_name %in% genesToMapToGO.base ]))
  
  genesToMapToGO.stress <- unlist(checkWhyNotEnriched[i, ]$Interacting_pairs_stress)
  df_stress = data.frame(table(Gomap_localization$GOnumber[Gomap_localization$sys_name %in% genesToMapToGO.stress ]))
  
  if(nrow(df_base) >0 & nrow(df_stress)>0){
    colnames(df_base) = c('GO', 'base_freq')
    colnames(df_stress) = c('GO', 'stress_freq')
    
    merged_df <- merge(df_base, df_stress, by.x='GO', by.y='GO', all.x=T, all.y=T)
    merged_df[is.na(merged_df)] = 0 
    print(merge(merged_df, GO_pathway_dict, by.x='GO', by.y='GOnumber', all.x=T, all.y=F))
    cat("\n")
    cat('-------------------------------------------------')
    cat("\n")
  }else{
    print('not included')
  }
    
}


### import GO map
GOmap <- read.delim('Data/GO/go_slim_mapping.tab', header = F)
colnames(GOmap)<- c('sys_name','stand_name','SGD_id', 'CFP','Pathway', 'GOnumber', 'type')
Gomap_localization <- GOmap[GOmap$CFP=='C', ]
table(Gomap_localization$GOnumber=='')

GO_pathway_dict <- unique(GOmap[,c('Pathway', 'GOnumber')])

## import Brown lab data
MMS.DataBrown <- read.csv('Data/Sources/localization_change_brown.csv')
MMS.DataBrown <- MMS.DataBrown[,! colnames(MMS.DataBrown) %in% c('X', 'X.1', 'X.2', 'X.3')]
head(MMS.DataBrown)

mmsRes <- read.csv('../Desktop/mmsRes.csv',stringsAsFactors = F)
head(mmsRes)

stress_base_ORFS_list <- readRDS('Results/listOf_localPartenrs_stress_base.rds')
MMS_base_ORFS <- stress_base_ORFS_list$MMS
tail(MMS_base_ORFS)
sum(MMS_base_ORFS$enriched_term != '')


stress_base_ORFS_localizPair <- readRDS('Results/listOf_stress_base_ORFS_localizPair.rds')
mms_base_ORFS_localizPair <- stress_base_ORFS_localizPair$MMS
tail(mms_base_ORFS_localizPair)

# how many of that list in generally included in ours
# how many of their changed localization are enriched ??
## investigate the differences

## 25% their proteins is included in our data , count= 63 -> 
##  only 1 is labeled as location change in our data
## we'll check the other 62 proteins in the following to see if they have actually 
### changed their location or not( maybe the fisher test is not working well)


sum(unique(MMS.DataBrown$Systematic.ORF) %in% 
      unique(MMS_base_ORFS$ORF))/length(unique(MMS.DataBrown$Systematic.ORF))

MMS.DataBrown_included <- MMS.DataBrown[MMS.DataBrown$Systematic.ORF %in% MMS_base_ORFS$ORF,]

### only 1 -> YJR056C
MMS_base_ORFS[MMS_base_ORFS$ORF %in% MMS.DataBrown_included$Systematic.ORF & MMS_base_ORFS$enriched_term  != '',]
mmsRes[mmsRes$ORF == 'YJR056C', ]
MMS.DataBrown_included[MMS.DataBrown_included$Systematic.ORF=='YJR056C',]

# GO:0005634	YJR056C
mms_base_ORFS_localizPair[which(MMS_base_ORFS$ORF == 'YJR056C' )]



#### these are the proteins that are changing their localization accourding to
####  the brown paper, but are not marked as changes by the fisher exact test
####  we attempt to check their frequency in both conditions to see if it had a noticible change or not

checkWhyNotEnriched <- MMS_base_ORFS[MMS_base_ORFS$ORF %in% MMS.DataBrown_included$Systematic.ORF 
                                     & MMS_base_ORFS$enriched_term  == '',]

sink("localization_comparison_BrownResults.txt")
sapply(1:nrow(checkWhyNotEnriched), function(i){
  get.Stress.Base.freq.table(i)
} )
sink()


###########
# show dk the ****** PTM results****** -> what is next ??
###########


## ToDO
# add the pre-calculated p-value to the table
# fix the Global frequency calculation
# calculate 2 other p-values for comparing the GO (or the protein??) with the global amount 
# check the code

# Annotate the proteins in the baseline condition by their neigbours -> find the localization -> pvalue by fisher exact test
##  -> Check if the GO annotation of proteins is consistent with this method of localization annotation
###  -> Check the HURI paper -> what was their approach for solving this 

# Compare with branda andrew lab results with HU with the mms results (same as the brown data)

## summerize the resuts and put them together -> make slides -> transfer the final results to DK




