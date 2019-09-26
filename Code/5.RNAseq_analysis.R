source('Code/Functions.R')
Initialize()

expression<- read.csv('Data/genes.fpkm_tracking_clean.csv')
expression$MMS_P_value[expression$MMS_P_value=='#DIV/0!'] <- 1
expression$Poor.Carbon_P_value[expression$Poor.Carbon_P_value=='#DIV/0!'] <- 1
expression$H2O2_P_value[expression$H2O2_P_value=='#DIV/0!'] <- 1
# write.csv(expression, 'Data/genes.fpkm_tracking_clean.csv')


FPKM_columns <- c("Baseline_01_FPKM", "Baseline_02_FPKM", "MMS_01_FPKM" , "MMS_02_FPKM", 
                  "PoorCarbon_01_FPKM" ,"PoorCarbon_02_FPKM" , "H2O2_01_FPKM" ,"H2O2_02_FPKM")
expression_FPKM <- expression[,colnames(expression) %in% FPKM_columns]
expression_FPKM_melt <- melt(expression_FPKM)
expression_FPKM_melt$variable <- as.character(expression_FPKM_melt$variable)
expression_FPKM_melt$variable <- gsub('_FPKM','',expression_FPKM_melt$variable)

p1=ggplot(expression_FPKM_melt, aes(x=value+0.1, color= variable))+
  geom_density()+theme_bw()+scale_x_log10()+xlab('FPKM+0.1(log)')+ylab('Sample')
p2=ggplot(expression_FPKM_melt, aes(x =variable, y=value+0.1, fill= variable))+
  geom_boxplot()+theme_bw()+scale_y_log10('FPKM+0.1(log)')+xlab('Sample')

head(expression_FPKM)
rownames(expression_FPKM) <- expression$tracking_id
colnames(expression_FPKM) <-  gsub('_FPKM','',colnames(expression_FPKM))
pheatmap(cor(expression_FPKM), main = 'mRNA-seq', color = inferno(10000),)

head(expression_FPKM)
clust_data <- t(expression_FPKM)
exp_PCA <- prcomp(clust_data, scale = F)
exp_PCA <- data.frame(subset(exp_PCA$x, select=c(PC1, PC2)))
exp_PCA$Group <- gsub('_.*', '', rownames(exp_PCA))
ggplot(exp_PCA, aes(x=PC1,y=PC2, color=Group))+geom_point(size=3)+
  theme_bw()+ggtitle('Principal Component Analysis- mRNAseq')

head(expression_FPKM)

####### Differential Gene Expression
rownames(expression_FPKM) <- expression$tracking_id
group <- gsub('_.*', '', colnames(expression_FPKM))
gr <- factor(group)
design <- model.matrix(~0+gr)
colnames(design) <- gsub('gr', '',colnames(design))

contr.matrix <- makeContrasts(
  BaseVsH2O2 = Baseline-H2O2, 
  BaseVsMMS = Baseline-MMS, 
  BaseVsPoorCarbon = Baseline-PoorCarbon, 
  H2O2vsMMS = H2O2-MMS,
  H2O2vsPoorCarbon = H2O2-PoorCarbon,
  MMMSvsPoorCarbon = MMS-PoorCarbon,
  levels = colnames(design))
contr.matrix


par(mfrow=c(1,2))
v <- voom(expression_FPKM, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")  
summary(decideTests(efit))
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

Base.vs.H2O2 = topTreat(tfit, coef=1, n=Inf)
Base.vs.MMS = topTreat(tfit, coef=2, n=Inf)
Base.vs.PoorCarbon = topTreat(tfit, coef=3, n=Inf)
H2O2.vs.MMS = topTreat(tfit, coef=4, n=Inf)
H2O2.vs.PoorCarbon = topTreat(tfit, coef=5, n=Inf)
MMS.vs.PoorCarbon = topTreat(tfit, coef=6, n=Inf)

## make color palette
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, 5))

#### Volcano Plots
p1=ggplot(Base.vs.H2O2, aes(x=logFC, y=-log10(P.Value)))+geom_point(aes(color=-log10(adj.P.Val)))+theme_bw()+ggtitle('Base.vs.H2O2')+sc
p2=ggplot(Base.vs.MMS, aes(x=logFC, y=-log10(P.Value)))+geom_point(aes(color=-log10(adj.P.Val)))+theme_bw()+ggtitle('Base.vs.MMS')+sc
p3=ggplot(Base.vs.PoorCarbon, aes(x=logFC, y=-log10(P.Value)))+geom_point(aes(color=-log10(adj.P.Val)))+theme_bw()+ggtitle('Base.vs.PoorCarbon')+sc
p4=ggplot(H2O2.vs.MMS, aes(x=logFC, y=-log10(P.Value)))+geom_point(aes(color=-log10(adj.P.Val)))+theme_bw()+ggtitle('H2O2.vs.MMS')+sc
p5=ggplot(H2O2.vs.PoorCarbon, aes(x=logFC, y=-log10(P.Value)))+geom_point(aes(color=-log10(adj.P.Val)))+theme_bw()+ggtitle('H2O2.vs.PoorCarbon')+sc
p6=ggplot(MMS.vs.PoorCarbon, aes(x=logFC, y=-log10(P.Value)))+geom_point(aes(color=-log10(adj.P.Val)))+theme_bw()+ggtitle('MMS.vs.PoorCarbon')+sc
pdf('Results/Plots/13.volcano_Plots.pdf', height = 9, width= 15)
grid.arrange(p1,p2,p3,p4,p5,p6, ncol=3, nrow=2)
dev.off()


listOf_tfitTab <- list(Base.vs.H2O2, Base.vs.MMS, Base.vs.PoorCarbon, H2O2.vs.MMS, H2O2.vs.PoorCarbon, MMS.vs.PoorCarbon)
names(listOf_tfitTab) <- c('Base.vs.H2O2', 'Base.vs.MMS', 'Base.vs.PoorCarbon', 'H2O2.vs.MMS', 'H2O2.vs.PoorCarbon', 'MMS.vs.PoorCarbon')
saveRDS(listOf_tfitTab, 'listOf_tfitTab_DE_analysis.rds')
listOf_DE <- lapply(listOf_tfitTab, function(x) subset(x, (logFC>2 | logFC< (-2)) & (P.Value<0.05) & (adj.P.Val<0.05)  ))

### ratio of DE genes/total genes
DE_ratio <- sapply(1:length(listOf_DE), function(i) 100*round(nrow(listOf_DE[[i]])/nrow(listOf_tfitTab[[i]]), 3))
DE_ratio.df <- data.frame(DE_ratio=DE_ratio, Samples= names(listOf_tfitTab) )

### Finding the gene name for the DE tracking IDs
listOf_DE_id <- lapply(listOf_DE, function(x) rownames(x))
listOf_DE_genes <- lapply(listOf_DE_id, function(x) expression$gene_short_name[expression$tracking_id %in% x ])
All_DE_Genes <- unlist(lapply(listOf_DE_genes, function(x) str_split(x, ',')))


ListOfHubs_comb <- readRDS('ListOfHubs_comb.rds')
ListOfHubs_comb2 <- lapply(ListOfHubs_comb, function(x) merge(x, GOterms, by.x='Gene', by.y='sys_name',all.x=T))
lapply(ListOfHubs_comb2, function(x) round(sum(x$stand_name %in% All_DE_Genes)/nrow(x),2))

ListOf_DE_Hubs <- lapply(ListOfHubs_comb2, function(x) x$stand_name[x$stand_name %in% All_DE_Genes])

vennPlot_DE_hubs <- venn.diagram(ListOf_DE_Hubs , NULL, fill=rainbow(4),
                            alpha=c(0.4,0.4,0.4,0.4), 
                            cex =2.2, cat.fontface=1, 
                            category.names=names(ListOf_DE_Hubs), 
                            main = 'Compare DE hubs between conditions')
grid.draw(vennPlot_DE_hubs)
