#!/usr/bin/env Rscript

source('Functions.R')
Initialize()

expression <- read.csv('genes.fpkm_tracking_clean.csv')
head(expression)

### Input data
GOterms <- read.delim('GO_terms_results.tsv', header = F)
colnames(GOterms) <- c('SGD_id', 'sys_name', 'Organism', 'stand_name', 'Description')

ListOfConditions <- readRDS('ListOfConditions_IS7.rds')
lapply(ListOfConditions, head)

AD_degrees <- lapply(ListOfConditions, .getAD_Degrees)
DB_degrees <- lapply(ListOfConditions, .getDB_Degrees)

AD_Genes <- lapply(AD_degrees, function(x) data.frame(Gene=x$Gene))
DB_Genes <- lapply(DB_degrees, function(x) data.frame(Gene=x$Gene))


############## Annotation issue
all_possible_genes <- data.frame(Genes=unique(c(as.character(unlist(AD_Genes)),as.character(unlist(DB_Genes)))))
all_possible_genes <- merge(all_possible_genes, GOterms, by.x='Genes', by.y='sys_name')
head(all_possible_genes)

listOfNames <- str_split(expression$gene_short_name, ',')

addColumns <- function(x){
  x$locus <- ''
  x$MMS_P_value <- -1        
  x$Poor.Carbon_P_value <- -1  
  x$H2O2_P_value <- -1 
  x$Baseline_01_FPKM <- -1     
  x$Baseline_02_FPKM <- -1  
  x$MMS_01_FPKM <- -1         
  x$MMS_02_FPKM <-  -1 
  x$PoorCarbon_01_FPKM <- -1 
  x$PoorCarbon_02_FPKM <- -1    
  x$H2O2_01_FPKM <- -1  
  x$H2O2_02_FPKM <- -1  
  return(x)
}

all_possible_genes <- addColumns(all_possible_genes)
head(all_possible_genes)

listOfUnmapped <- c()
for (i in 1:nrow(all_possible_genes)){
  print(i)
  flag=FALSE
  for(j in 1:nrow(expression)){
    if( (all_possible_genes[i,]$Genes %in% listOfNames[[j]] & all_possible_genes[i,]$Genes != '') | 
       (all_possible_genes[i,]$stand_name %in% listOfNames[[j]] & all_possible_genes[i,]$stand_name != '')){
      flag=TRUE
      all_possible_genes[i,]$locus <- expression[j,]$locus
      all_possible_genes[i,]$MMS_P_value <- expression[j,]$MMS_P_value       
      all_possible_genes[i,]$Poor.Carbon_P_value <- expression[j,]$Poor.Carbon_P_value 
      all_possible_genes[i,]$H2O2_P_value <- expression[j,]$H2O2_P_value 
      all_possible_genes[i,]$Baseline_01_FPKM <- expression[j,]$Baseline_01_FPKM     
      all_possible_genes[i,]$Baseline_02_FPKM <- expression[j,]$Baseline_02_FPKM  
      all_possible_genes[i,]$MMS_01_FPKM <- expression[j,]$MMS_01_FPKM         
      all_possible_genes[i,]$MMS_02_FPKM <-  expression[j,]$MMS_02_FPKM 
      all_possible_genes[i,]$PoorCarbon_01_FPKM <- expression[j,]$PoorCarbon_01_FPKM
      all_possible_genes[i,]$PoorCarbon_02_FPKM <- expression[j,]$PoorCarbon_02_FPKM    
      all_possible_genes[i,]$H2O2_01_FPKM <- expression[j,]$H2O2_01_FPKM 
      all_possible_genes[i,]$H2O2_02_FPKM <- expression[j,]$H2O2_02_FPKM 
    }
  }
  if(!flag) listOfUnmapped <- c(listOfUnmapped, i)
}

saveRDS(all_possible_genes, 'all_possible_genes_exp.rds')
saveRDS(listOfUnmapped, 'listOfUnmapped.rds')

