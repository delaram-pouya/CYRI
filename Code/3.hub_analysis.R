## import needed functions and libraries
source('Code/Functions.R')
Initialize()
## given a list of hubs for each consition, compares 
##  how many of these hubs are shared between different conditions


### Input data
ListOfConditions <- readRDS('Data/ListOfConditions_IS7.rds')
lapply(ListOfConditions, head)

AD_degrees <- lapply(ListOfConditions, .getAD_Degrees)
DB_degrees <- lapply(ListOfConditions, .getDB_Degrees)

AD_genes <- lapply(AD_degrees, function(x) unique(x$Gene))
DB_genes <- lapply(DB_degrees, function(x) unique(x$Gene))
sapply(1:length(DB_genes), function(i)sum(AD_genes[[i]] %in% DB_genes[[i]]))
lapply(DB_genes, length)


### compare AD degree distributions
AD_degrees_df <- sapply(1:length(AD_degrees), function(i) 
  data.frame(degree=AD_degrees[[i]]$Freq, data=names(ListOfConditions)[i]), simplify = F)
AD_degrees_df <- do.call(rbind, AD_degrees_df)

p_AD_degree <- ggplot(AD_degrees_df ,aes(x=degree, color= data))+geom_density()+
  theme_bw()+ggtitle('AD degree distribution')+scale_x_log10()


### comapre DB degree distributions
DB_degrees_df <- sapply(1:length(DB_degrees), function(i) 
  data.frame(degree=DB_degrees[[i]]$Freq, data=names(ListOfConditions)[i]), simplify = F)
DB_degrees_df <- do.call(rbind, DB_degrees_df)

p_DB_degree <- ggplot(DB_degrees_df ,aes(x=degree, color= data))+geom_density()+
  theme_bw()+ggtitle('DB degree distribution')+scale_x_log10()

grid.arrange(p_AD_degree, p_DB_degree, nrow=1, ncol= 2)

### how to combine AD and DB degrees????

### Find the hubs
###  based on 95% threshold

### AD 
hub_thresholds_AD <- unlist(lapply(AD_degrees, function(x) quantile(x$Freq, 0.95)))
ListOfHubs_AD <- sapply(1:length(AD_degrees), function(i) subset(AD_degrees[[i]], Freq>hub_thresholds_AD[i]), simplify = F)
names(ListOfHubs_AD) <- names(ListOfConditions)
sapply(1:length(ListOfHubs_AD), function(i) paste0(names(ListOfConditions)[i],': ',nrow(ListOfHubs_AD[[i]]),'/',nrow(AD_degrees[[i]])))
saveRDS(ListOfHubs_AD , 'Results/ListOfHubs_AD.rds')


### DB
hub_thresholds_DB <- unlist(lapply(DB_degrees, function(x) quantile(x$Freq, 0.95)))
ListOfHubs_DB <- sapply(1:length(DB_degrees), function(i) subset(DB_degrees[[i]], Freq>hub_thresholds_DB[i]), simplify = F)
names(ListOfHubs_DB) <- names(ListOfConditions)
sapply(1:length(ListOfHubs_DB), function(i) paste0(names(ListOfConditions)[i],': ',nrow(ListOfHubs_DB[[i]]),'/',nrow(DB_degrees[[i]])))
saveRDS(ListOfHubs_DB , 'Results/ListOfHubs_DB.rds')


compare_hubs_AD <- getHubComparisonMatrix(ListOfHubs_AD)
p_comp_hubs_AD <- ggplot(compare_hubs_AD, aes(x=sum))+geom_bar(color='black', fill='pink2', alpha=0.7)+
  theme_bw()+xlab('#shared-hubs')+ggtitle('Compare AD hubs in different conditions')
p_comp_hubs_AD
round(table(compare_hubs_AD$sum)*100/sum(table(compare_hubs_AD$sum)), 2)


compare_hubs_DB <- getHubComparisonMatrix(ListOfHubs_DB)
p_comp_hubs_DB <- ggplot(compare_hubs_DB, aes(x=sum))+geom_bar(color='black', fill='green4', alpha=0.7)+
  theme_bw()+xlab('#shared-hubs')+ggtitle('Compare DB hubs in different conditions')
p_comp_hubs_DB
round(table(compare_hubs_DB$sum)*100/sum(table(compare_hubs_DB$sum)), 2)

grid.arrange(p_comp_hubs_AD, p_comp_hubs_DB, nrow=1, ncol= 2)

ListOfHubs_AD_gene <- lapply(ListOfHubs_AD, function(x) as.character(x$Gene))
vennPlot_AD <- venn.diagram(ListOfHubs_AD_gene , NULL, fill=c("cyan1", "deeppink4",'gold','darkblue'),
             alpha=c(0.4,0.4,0.4,0.4), 
             cex =2.2, cat.fontface=1, category.names=names(test), main = 'Hubs-AD')
grid.draw(vennPlot_AD)



#### Venn Diagrams
ListOfHubs_AD_gene <- lapply(ListOfHubs_AD, function(x) as.character(x$Gene))
vennPlot_AD <- venn.diagram(ListOfHubs_AD_gene , NULL, fill=c("cyan1", "deeppink4",'gold','darkblue'),
                            alpha=c(0.4,0.4,0.4,0.4), 
                            cex =2.2, cat.fontface=1, category.names=names(test), main = 'Hubs-AD')
grid.draw(vennPlot_AD)


ListOfHubs_DB_gene <- lapply(ListOfHubs_DB, function(x) as.character(x$Gene))
vennPlot_DB <- venn.diagram(ListOfHubs_DB_gene , NULL, fill=c("cyan1", "deeppink4",'gold','darkblue'),
                            alpha=c(0.4,0.4,0.4,0.4), 
                            cex =2.2, cat.fontface=1, category.names=names(test), main = 'Hubs-DB')
grid.draw(vennPlot_DB)



####### combine the AD and DB information
### first get the unique PPIs so that the design bias is removed

temp <- lapply(ListOfConditions, function(x) unique(x$PPI))
temp <- lapply(temp, function(x) str_split_fixed(x, '_', 2))
comb_degrees <- lapply(temp, function(x) data.frame(table(x)))
comb_degrees <- lapply(comb_degrees, function(x) {x=x[order(x$Freq, decreasing = T),]; colnames(x)[1]='Gene';x})
saveRDS(comb_degrees, 'comb_degrees.rds')
lapply(comb_degrees, head)

hub_thresholds_comb <- unlist(lapply(comb_degrees, function(x) quantile(x$Freq, 0.95)))
ListOfHubs_comb <- sapply(1:length(comb_degrees), function(i) subset(comb_degrees[[i]], Freq>hub_thresholds_comb[i]), simplify = F)
names(ListOfHubs_comb) <- names(ListOfConditions)
lapply(ListOfHubs_comb, print)

compare_hubs_comb <- getHubComparisonMatrix(ListOfHubs_comb)
p_comp_hubs_comb <- ggplot(compare_hubs_comb, aes(x=sum))+geom_bar(color='black', fill='cyan3', alpha=0.7)+
  theme_bw()+xlab('#shared-hubs')+ggtitle('Compare combined hubs in different conditions')
p_comp_hubs_comb

ListOfHubs_comb_gene <- lapply(ListOfHubs_comb, function(x) as.character(x$Gene))
vennPlot_comb <- venn.diagram(ListOfHubs_comb_gene , NULL, fill=c("cyan1", "deeppink4",'gold','darkblue'),
                            alpha=c(0.4,0.4,0.4,0.4), 
                            cex =2.2, cat.fontface=1, category.names=names(ListOfHubs_comb_gene), main = 'Hubs-combined')
grid.draw(vennPlot_comb)
# saveRDS(ListOfHubs_comb, 'ListOfHubs_comb.rds')


###############  make a dataframe of AD and DB protein degrees in each condition 
AD_degree_Matrix <- getDegreeMatrix(AD_degrees)
DB_degree_Matrix <- getDegreeMatrix(DB_degrees)
head(AD_degree_Matrix)
head(DB_degree_Matrix)
AD_degree_Matrix <- readRDS('Results/AD_degree_Matrix.rds')
DB_degree_Matrix <- readRDS('Results/DB_degree_Matrix.rds')

total_degree_matrix <- merge(AD_degree_Matrix, DB_degree_Matrix, by.x='Gene', by.y= 'Gene', all=T)
### many of the proteins are not present at all in the different conditions ?????
head(total_degree_matrix)


########### Manual check of interesting conditions

###### AD 
## Proteins which were not a hub in the baseline data but are a hub in the changed condtion
subset(compare_hubs_AD, Baseline==0)
input_df <- subset(compare_hubs_AD, Baseline==0 & sum==3)
resMAT_1 <- get_attributes(subset(compare_hubs_AD, Baseline==0 & sum>=2), AD_degree_Matrix) ## interesting!!! HAT1 :)
write.csv(resMAT_1, 'Results/Result_Tables/1.stress_higher_than_baseline_AD.csv', quote = F, row.names = F)

## proteins which degree has decreased in specific conditions
resMAT_2 <- get_attributes(subset(compare_hubs_AD, MMS==0 & sum==3), AD_degree_Matrix)
write.csv(resMAT_2, 'Results/Result_Tables/2.MMS_lower_than_all_AD.csv', quote = F, row.names = F)
resMAT_3 <- get_attributes(subset(compare_hubs_AD, PoorCarbon==0 & sum==3), AD_degree_Matrix)
write.csv(resMAT_3, 'Results/Result_Tables/3.PoorCarbon_lower_than_all_AD.csv', quote = F, row.names = F)

#### proteins which degree has increased in specific conditions
resMAT_4 <- get_attributes(subset(compare_hubs_AD, MMS==1 & sum==1), AD_degree_Matrix)
write.csv(resMAT_4, 'Results/Result_Tables/4.MMS_higher_than_all_AD.csv', quote = F, row.names = F)

resMAT_5 <- get_attributes(subset(compare_hubs_AD, H2O2==1 & sum==1), AD_degree_Matrix)
write.csv(resMAT_5, 'Results/Result_Tables/5.H2O2_higher_than_all_AD.csv', quote = F, row.names = F)

resMAT_6 <- get_attributes(subset(compare_hubs_AD, PoorCarbon==1 & sum==1), AD_degree_Matrix)
write.csv(resMAT_6, 'Results/Result_Tables/6.PoorCarbon_higher_than_all_AD.csv', quote = F, row.names = F)

### proteins which degree has decreased in all conditions
resMAT_7 <- get_attributes(subset(compare_hubs_AD, Baseline==1 & sum==1), AD_degree_Matrix)
write.csv(resMAT_7, 'Results/Result_Tables/7.baseline_higher_than_stress_AD.csv', quote = F, row.names = F)



############# DB

subset(compare_hubs_DB, Baseline==0)
input_df <- subset(compare_hubs_DB, Baseline==0 & sum==3)
resMAT_1 <- get_attributes(subset(compare_hubs_DB, Baseline==0 & sum>=2), DB_degree_Matrix) 
write.csv(resMAT_1, 'Results/Result_Tables/1.stress_higher_than_baseline_DB.csv', quote = F, row.names = F)

## proteins which degree has decreased in specific conditions
resMAT_2 <- get_attributes(subset(compare_hubs_DB, MMS==0 & sum==3), DB_degree_Matrix)
write.csv(resMAT_2, 'Results/Result_Tables/2.MMS_lower_than_all_DB.csv', quote = F, row.names = F)
resMAT_3 <- get_attributes(subset(compare_hubs_DB, PoorCarbon==0 & sum==3), DB_degree_Matrix)
write.csv(resMAT_3, 'Results/Result_Tables/3.PoorCarbon_lower_than_all_DB.csv', quote = F, row.names = F)

#### proteins which degree has increased in specific conditions
resMAT_4 <- get_attributes(subset(compare_hubs_DB, MMS==1 & sum==1), DB_degree_Matrix)
write.csv(resMAT_4, 'Results/Result_Tables/4.MMS_higher_than_all_DB.csv', quote = F, row.names = F)

resMAT_5 <- get_attributes(subset(compare_hubs_DB, H2O2==1 & sum==1), DB_degree_Matrix)
write.csv(resMAT_5, 'Results/Result_Tables/5.H2O2_higher_than_all_DB.csv', quote = F, row.names = F)

resMAT_6 <- get_attributes(subset(compare_hubs_DB, PoorCarbon==1 & sum==1), DB_degree_Matrix)
write.csv(resMAT_6, 'Results/Result_Tables/6.PoorCarbon_higher_than_all_DB.csv', quote = F, row.names = F)

### proteins which degree has decreased in all conditions
resMAT_7 <- get_attributes(subset(compare_hubs_DB, Baseline==1 & sum==1), DB_degree_Matrix)
write.csv(resMAT_7, 'Results/Result_Tables/7.baseline_higher_than_stress_DB.csv', quote = F, row.names = F)


