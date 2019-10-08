source('Code/Functions.R')
Initialize()
ListOfConditions <- readRDS('Data/ListOfConditions_IS7.rds')
lapply(ListOfConditions, head)
#####################

## Functions 
get.AD_DB.list <- function(AD_degrees.condition, DB_degrees.condition){
  AD_DB <- list(as.character(AD_degrees.condition$Gene), as.character(DB_degrees.condition$Gene))
  names(AD_DB) <- c('AD', 'DB')
  return(AD_DB)
}


get.Venn.plot <- function(listToDraw, Title){
  venn.diagram( listToDraw , NULL, fill=rainbow(2),
                alpha=c(0.4,0.4), 
                cex =2.3, cat.fontface=1, 
                category.names=names(listToDraw), main = Title , main.cex = 2)
}



## check degree correlation between the shared proteins in AD and DB score
## duplicates (ppi) are not removed
AD_degrees <- lapply(ListOfConditions, .getAD_Degrees)
DB_degrees <- lapply(ListOfConditions, .getDB_Degrees)
lapply(AD_degrees, head)
lapply(DB_degrees, head)

ListOfHubs_DB <- readRDS( 'Results/ListOfHubs_DB.rds')
ListOfHubs_AD <- readRDS( 'Results/ListOfHubs_AD.rds')
lapply(ListOfHubs_AD, head)
lapply(ListOfHubs_DB, head)

## how many(percentage) of the PPI pairs are deplicates ?
lapply(ListOfConditions, 
       function(x){
         df = data.frame(table(x$PPI))
         round(table(df$Freq)/sum(table(df$Freq)), 2)})


#### Compare the AD/DB for each condition
baseline.AD_DB <- get.AD_DB.list(AD_degrees$Baseline, DB_degrees$Baseline)
grid.draw(get.Venn.plot(baseline.AD_DB, 'Baseline AD/DB'))

h2o2.AD_DB <- get.AD_DB.list(AD_degrees$H2O2, DB_degrees$H2O2)
grid.draw(get.Venn.plot(h2o2.AD_DB, 'H2O2 AD/DB'))

mms.AD_DB <- get.AD_DB.list(AD_degrees$MMS, DB_degrees$MMS)
grid.draw(get.Venn.plot(mms.AD_DB, 'MMS AD/DB'))

poorcarbon.AD_DB <- get.AD_DB.list(AD_degrees$PoorCarbon, DB_degrees$PoorCarbon)
grid.draw(get.Venn.plot(poorcarbon.AD_DB, 'Poor Carbon AD/DB'))





  
  
  
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

