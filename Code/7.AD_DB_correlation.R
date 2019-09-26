source('Code/Functions.R')
Initialize()
ListOfConditions <- readRDS('Data/ListOfConditions_IS7.rds')
lapply(ListOfConditions, head)
#####################
## check degree correlation between the shared proteins in AD and DB score

AD_degrees <- lapply(ListOfConditions, .getAD_Degrees)
DB_degrees <- lapply(ListOfConditions, .getDB_Degrees)
lapply(AD_degrees, head)
lapply(DB_degrees, head)

ListOfHubs_DB <- readRDS( 'Results/ListOfHubs_DB.rds')
ListOfHubs_AD <- readRDS( 'Results/ListOfHubs_AD.rds')
lapply(ListOfHubs_AD, head)
lapply(ListOfHubs_DB, head)


AD_DB_shared <- sapply(1:length(AD_degrees), function(i) {
  x = merge(AD_degrees[[i]], DB_degrees[[i]], by.x='Gene', by.y='Gene',all.x=F, all.y=F)
  colnames(x) <- c('Gene', 'AD_degree','DB_degree');x
  },simplify = F)

names(AD_DB_shared) <- names(ListOfConditions)

#### split the list
Baseline_shared <- AD_DB_shared$Baseline
H2O2_shared <- AD_DB_shared$H2O2
MMS_shared <- AD_DB_shared$MMS
PoorCarbon_shared <- AD_DB_shared$PoorCarbon
head(PoorCarbon_shared)

## Baseline
Baseline_shared$AD_hub <- ifelse(Baseline_shared$Gene %in% ListOfHubs_AD$Baseline$Gene, 1 , 0)
Baseline_shared$DB_hub <- ifelse(Baseline_shared$Gene %in% ListOfHubs_DB$Baseline$Gene, 1 , 0)
Baseline_shared$HUB <- ifelse(Baseline_shared$AD_hub==1 & Baseline_shared$DB_hub==1, 'Both', 
                              ifelse(Baseline_shared$AD_hub==1 & Baseline_shared$DB_hub==0, 'AD',
                                     ifelse(Baseline_shared$AD_hub==0 & Baseline_shared$DB_hub==1, 'DB','Non')))
  
## H2O2 
H2O2_shared$AD_hub <- ifelse(H2O2_shared$Gene %in% ListOfHubs_AD$H2O2$Gene, 1 , 0)
H2O2_shared$DB_hub <- ifelse(H2O2_shared$Gene %in% ListOfHubs_DB$H2O2$Gene, 1 , 0)
H2O2_shared$HUB <- ifelse(H2O2_shared$AD_hub==1 & H2O2_shared$DB_hub==1, 'Both', 
                              ifelse(H2O2_shared$AD_hub==1 & H2O2_shared$DB_hub==0, 'AD',
                                     ifelse(H2O2_shared$AD_hub==0 & H2O2_shared$DB_hub==1, 'DB','Non')))

## MMS
MMS_shared$AD_hub <- ifelse(MMS_shared$Gene %in% ListOfHubs_AD$MMS$Gene, 1 , 0)
MMS_shared$DB_hub <- ifelse(MMS_shared$Gene %in% ListOfHubs_DB$MMS$Gene, 1 , 0)
MMS_shared$HUB <- ifelse(MMS_shared$AD_hub==1 & MMS_shared$DB_hub==1, 'Both', 
                              ifelse(MMS_shared$AD_hub==1 & MMS_shared$DB_hub==0, 'AD',
                                     ifelse(MMS_shared$AD_hub==0 & MMS_shared$DB_hub==1, 'DB','Non')))

## PoorCarbon
PoorCarbon_shared$AD_hub <- ifelse(PoorCarbon_shared$Gene %in% ListOfHubs_AD$PoorCarbon$Gene, 1 , 0)
PoorCarbon_shared$DB_hub <- ifelse(PoorCarbon_shared$Gene %in% ListOfHubs_DB$PoorCarbon$Gene, 1 , 0)
PoorCarbon_shared$HUB <- ifelse(PoorCarbon_shared$AD_hub==1 & PoorCarbon_shared$DB_hub==1, 'Both', 
                              ifelse(PoorCarbon_shared$AD_hub==1 & PoorCarbon_shared$DB_hub==0, 'AD',
                                     ifelse(PoorCarbon_shared$AD_hub==0 & PoorCarbon_shared$DB_hub==1, 'DB','Non')))

pdf('AD_DB_hub_scatter_plot.pdf', height = 5, width = 5)
ggplot(Baseline_shared, aes(log2(AD_degree), log2(DB_degree), color=HUB))+geom_point()+theme_bw()+ggtitle('Baseline')
ggplot(H2O2_shared, aes(log2(AD_degree), log2(DB_degree), color=HUB))+geom_point()+theme_bw()+ggtitle('H2O2')
ggplot(MMS_shared, aes(log2(AD_degree), log2(DB_degree), color=HUB))+geom_point()+theme_bw()+ggtitle('MMS')
ggplot(PoorCarbon_shared, aes(log2(AD_degree), log2(DB_degree), color=HUB))+geom_point()+theme_bw()+ggtitle('PoorCarbon')
dev.off()


#pdf('AD_DB_degree_cor_log2.pdf')
sapply(1:length(AD_DB_shared), function(i){
  ggplot(AD_DB_shared[[i]], aes(x=log2(AD_degree), y=log2(DB_degree)))+
    geom_point(size=2)+theme_bw()+ggtitle(names(AD_DB_shared)[i])
}, simplify = F)
dev.off()


###### Pearson correlation
pearson <- sapply(1:length(AD_DB_shared), function(i){
  cor.test(AD_DB_shared[[i]]$AD_degree, AD_DB_shared[[i]]$DB_degree, method = 'pearson')
}, simplify = F)
names(pearson) <- names(ListOfConditions)
lapply(pearson, function(x) x$estimate)
lapply(pearson, function(x) x$p.value)

###### Spearman Correlation
spearman <- sapply(1:length(AD_DB_shared), function(i){
  cor.test(AD_DB_shared[[i]]$AD_degree, AD_DB_shared[[i]]$DB_degree, method = 'spearman',exact=FALSE)
}, simplify = F)
names(spearman) <- names(ListOfConditions)
lapply(spearman, function(x) x$estimate)
lapply(spearman, function(x) x$p.value)


lapply(AD_DB_shared, nrow)
lapply(AD_degrees, nrow)
lapply(DB_degrees, nrow)

