# Check the 2 scenarios DK mentioned and think about a way to detect them
# Endogenous gene higher expression -> competition -> lower signal
# A third gene over-expression -> high signal  

#### Need to run/check the exp mapping again -> some genes were not considered
##### considering that the matrix is complete -> how to address the Qs that DK mentioned before ?? 

source('Code/Functions.R')
Initialize()
library(reshape2)
library("gplots")
library(corrplot)

## Functions
getExpStatus <- function( exp_difference_with_base , setOf_exp_stress_base){
  ifelse(exp_difference_with_base < quantile(setOf_exp_stress_base, 0.25) | 
           exp_difference_with_base == quantile(setOf_exp_stress_base, 0.25), 'eLESS', 
         ifelse(exp_difference_with_base < quantile(setOf_exp_stress_base, 0.75), 'eEQUAL', 'eMORE'))
}

getDegStatus <- function( deg_difference_with_base ,setOf_deg_stress_base){
  ifelse(deg_difference_with_base < quantile(setOf_deg_stress_base, 0.25) | 
           deg_difference_with_base == quantile(setOf_deg_stress_base, 0.25), 'dLESS', 
         ifelse(deg_difference_with_base < quantile(setOf_deg_stress_base, 0.75), 'dEQUAL', 'dMORE'))
}


getEndogensExpDegreeMatrix <- function(ppi_exp_stress){
  matrix <- data.frame(exp=c(rep('eLESS', 3), rep('eEQUAL', 3), rep('eMORE', 3)), 
                       degree=rep(c('dLESS', 'dEQUAL', 'dMORE'), 3))
  
  matrix$count <- c(nrow(subset(ppi_exp_stress, exp_status == 'eLESS' & deg_status == 'dLESS')), 
                    nrow(subset(ppi_exp_stress, exp_status == 'eLESS' & deg_status == 'dEQUAL')),
                    nrow(subset(ppi_exp_stress, exp_status == 'eLESS' & deg_status == 'dMORE')),
                    
                    nrow(subset(ppi_exp_stress, exp_status == 'eEQUAL' & deg_status == 'dLESS')), 
                    nrow(subset(ppi_exp_stress, exp_status == 'eEQUAL' & deg_status == 'dEQUAL')),
                    nrow(subset(ppi_exp_stress, exp_status == 'eEQUAL' & deg_status == 'dMORE')),
                    
                    nrow(subset(ppi_exp_stress, exp_status == 'eMORE' & deg_status == 'dLESS')), 
                    nrow(subset(ppi_exp_stress, exp_status == 'eMORE' & deg_status == 'dEQUAL')),
                    nrow(subset(ppi_exp_stress, exp_status == 'eMORE' & deg_status == 'dMORE')))
  return(matrix)
}


expression<- read.csv('Data/genes.fpkm_tracking_clean.csv')
all_possible_genes_filled <- readRDS('mapGeneExp/all_possible_genes_exp.rds')
all_possible_genes_filled[all_possible_genes_filled==-1] <- NA
unmaped_indices <- readRDS('mapGeneExp/listOfUnmapped.rds')

### Input data
ListOfConditions <- readRDS('Data/ListOfConditions_IS7.rds')
lapply(ListOfConditions, head)

### you get very different results by counting or not counting the duplicates
combined_degrees <- lapply(ListOfConditions, .getComb_Degree)
lapply(combined_degrees, head)
combined_degree_matrix <- getDegreeMatrix(combined_degrees)

all_possible_genes_filled = merge(all_possible_genes_filled, combined_degree_matrix, by.x='Genes', by.y='Gene', all.x=T)
head(all_possible_genes_filled)

### Scenerio number 1 
### endogenous higher expression in the new condition -> 
#### lead to competetion with our product -> lower signal

#### mean of the expression is considered
ppi_exp_merged <- data.frame( standard_name= all_possible_genes_filled$stand_name,
                              systematic_name= all_possible_genes_filled$Genes,
                              
                              exp_Baseline=rowMeans(cbind(all_possible_genes_filled$Baseline_01_FPKM, all_possible_genes_filled$Baseline_02_FPKM), na.rm=TRUE), 
                              exp_MMS=rowMeans(cbind(all_possible_genes_filled$MMS_01_FPKM, all_possible_genes_filled$MMS_02_FPKM), na.rm=TRUE), 
                              exp_PoorCarbon=rowMeans(cbind(all_possible_genes_filled$PoorCarbon_01_FPKM, all_possible_genes_filled$PoorCarbon_02_FPKM), na.rm=TRUE),
                              exp_H2O2=rowMeans(cbind(all_possible_genes_filled$H2O2_01_FPKM, all_possible_genes_filled$H2O2_02_FPKM), na.rm=TRUE), 
                              
                              deg_Baseline = all_possible_genes_filled$Freq_Baseline, 
                              deg_H2O2 = all_possible_genes_filled$Freq_H2O2, 
                              deg_MMS = all_possible_genes_filled$Freq_MMS, 
                              deg_PoorCarbon = all_possible_genes_filled$Freq_PoorCarbon)

### is subtraction a good idea?? maybe getting the log2(devide) would work better ????
ppi_exp_merged <- ppi_exp_merged[complete.cases(ppi_exp_merged), ]
dim(ppi_exp_merged)

### calculate the expression difference between Stress condition and Baseline
ppi_exp_merged$Exp_mms_base= log2((ppi_exp_merged$exp_MMS+1e-4)/(ppi_exp_merged$exp_Baseline+1e-4))
ppi_exp_merged$Exp_poorcarbon_base = log2((ppi_exp_merged$exp_PoorCarbon+1e-4) / (ppi_exp_merged$exp_Baseline+1e-4))
ppi_exp_merged$Exp_h2o2_base = log2((ppi_exp_merged$exp_H2O2+1e-4) / (ppi_exp_merged$exp_Baseline+1e-4))

summary(ppi_exp_merged$Exp_poorcarbon_base)

### calculate the degree difference between the stress and Baseline
ppi_exp_merged$Deg_mms_base <-  ppi_exp_merged$deg_MMS - ppi_exp_merged$deg_Baseline
ppi_exp_merged$Deg_poorcarbon_base <-  ppi_exp_merged$deg_PoorCarbon - ppi_exp_merged$deg_Baseline
ppi_exp_merged$Deg_h2o2_base <-  ppi_exp_merged$deg_H2O2 - ppi_exp_merged$deg_Baseline

ppi_exp_merged <- readRDS('ppi_exp_merged.rds')


##### FIND the appropiate threshold here :
head(ppi_exp_merged)
setOf_exp_stress_base <- c(ppi_exp_merged$Exp_mms_base, ppi_exp_merged$Exp_h2o2_base, ppi_exp_merged$Exp_poorcarbon_base)
setOf_deg_stress_base <- c(ppi_exp_merged$Deg_mms_base, ppi_exp_merged$Deg_poorcarbon_base, ppi_exp_merged$Deg_h2o2_base)


########### MMS
ppi_exp_mms <- ppi_exp_merged[,colnames(ppi_exp_merged) %in% c('standard_name', 'Exp_mms_base', 'Deg_mms_base') ]
ppi_exp_mms$exp_status <- getExpStatus(ppi_exp_mms$Exp_mms_base, setOf_exp_stress_base)
ppi_exp_mms$deg_status <- getDegStatus(ppi_exp_mms$Deg_mms_base, setOf_deg_stress_base)
head(ppi_exp_mms)

## to check the specific examples:
subset(ppi_exp_mms, exp_status=='eMORE' &  deg_status=='dLESS')
subset(ppi_exp_mms, Exp_mms_base>0.4 &  Deg_mms_base<(-2))


#### POORCARBON is not relaible since it has high expression -> needs normalization
########### PoorCarbon
ppi_exp_poorcarbon <- ppi_exp_merged[,colnames(ppi_exp_merged) %in% c('standard_name','Exp_poorcarbon_base', 'Deg_poorcarbon_base') ]
ppi_exp_poorcarbon$exp_status <- getExpStatus(ppi_exp_poorcarbon$Exp_poorcarbon_base, setOf_exp_stress_base)
ppi_exp_poorcarbon$deg_status <- getDegStatus(ppi_exp_poorcarbon$Deg_poorcarbon_base, setOf_deg_stress_base)
head(ppi_exp_poorcarbon)

## to check the specific examples:
subset(ppi_exp_poorcarbon, exp_status=='eMORE' &  deg_status=='dLESS')

########### H2O2
ppi_exp_h2o2 <- ppi_exp_merged[,colnames(ppi_exp_merged) %in% c('standard_name', 'Exp_h2o2_base', 'Deg_h2o2_base') ]
ppi_exp_h2o2$exp_status <- getExpStatus(ppi_exp_h2o2$Exp_h2o2_base, setOf_exp_stress_base)
ppi_exp_h2o2$deg_status <- getDegStatus(ppi_exp_h2o2$Deg_h2o2_base, setOf_deg_stress_base)
head(ppi_exp_h2o2)

## to check the specific examples:
subset(ppi_exp_h2o2, exp_status=='eMORE' &  deg_status=='dLESS')




#####################
### Endogenous Gene Expression & Degree relationship
matrix_h2o2 <- getEndogensExpDegreeMatrix(ppi_exp_h2o2)
matrix_mms <- getEndogensExpDegreeMatrix(ppi_exp_mms)
matrix_poorcarbon <- getEndogensExpDegreeMatrix(ppi_exp_poorcarbon)

exp_deg_table_h2o2 <- acast(matrix_h2o2,exp~degree, value.var='count' )
exp_deg_table_mms <- acast(matrix_mms,exp~degree, value.var='count' )
exp_deg_table_poorcarbon <- acast(matrix_poorcarbon,exp~degree, value.var='count' )


pdf('ballonPlots.pdf')
balloonplot(t(as.table(as.matrix(exp_deg_table_h2o2))), 
            main ="H2O2 ", xlab ="Degree", ylab="Expression",
            label = TRUE, show.margins = FALSE, label.size=0.6, text.size=0.6)

balloonplot(t(as.table(as.matrix(exp_deg_table_mms))), 
            main ="MMS", xlab ="Degree", ylab="Expression",
            label = TRUE, show.margins = FALSE, label.size=0.6, text.size=0.6)

balloonplot(t(as.table(as.matrix(exp_deg_table_poorcarbon))), 
            main ="Poor-Carbon", xlab ="Degree", ylab="Expression",
            label = TRUE, show.margins = FALSE, label.size=0.6, text.size=0.6)
dev.off()
# top 5 



#############################
# http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r
chisq_h2o2 <- chisq.test(exp_deg_table_h2o2)
chisq_mms <- chisq.test(exp_deg_table_mms)
chisq_poorcarbon <- chisq.test(exp_deg_table_poorcarbon)
chisq$observed
chisq$expected
round(chisq$residuals, 3)
corrplot(chisq$residuals, is.cor = FALSE)
# Contibution in percentage (%)
contrib <- 100*chisq$residuals^2/chisq$statistic
round(contrib, 3)
# Visualize the contribution
corrplot(contrib, is.cor = FALSE)


################# Case-Study
head(ppi_exp_merged)
summary(setOf_exp_stress_base)
summary(setOf_deg_stress_base)

caseStudy_mms <- subset(ppi_exp_merged, Deg_mms_base< (-1) &  Exp_mms_base> 0.5)
caseStudy_mms <- subset(caseStudy_mms, select=c(standard_name ,systematic_name, exp_Baseline, exp_MMS, deg_Baseline, deg_MMS))
write.csv(caseStudy_mms, '../Desktop/caseStudy_mms.csv')

caseStudy_h2o2 <-subset(ppi_exp_merged, Deg_h2o2_base< (-1) &  Exp_h2o2_base> 0.4)
caseStudy_h2o2 <- subset(caseStudy_h2o2, select=c(standard_name ,systematic_name, exp_Baseline, exp_H2O2, deg_Baseline, deg_H2O2))
write.csv(caseStudy_h2o2, '../Desktop/caseStudy_h2o2.csv')

caseStudy_poorcarbon <-subset(ppi_exp_merged, Deg_poorcarbon_base< (-3) &  Exp_poorcarbon_base > 2)
caseStudy_poorcarbon <- subset(caseStudy_poorcarbon, select=c(standard_name ,systematic_name, exp_Baseline, exp_PoorCarbon, deg_Baseline, deg_PoorCarbon))
write.csv(caseStudy_poorcarbon, '../Desktop/caseStudy_poorcarbon.csv')


