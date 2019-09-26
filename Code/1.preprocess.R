## CYRI data analysis
## import needed functions and libraries
source('Code/Functions.R')
Initialize()

## import data
PATH <- 'Data/IS7'
files <- list.files(PATH, full.names = T, include.dirs = T)
data_S7 <- lapply(files, read.csv, header= F)
names(data_S7) <- str_split_fixed(list.files(PATH), '_', 2)[,1]
lapply(data_S7, head)

###  Cleaning the data
ADhub <- .cleanHub(data_S7[['ADhub']])
DBhub <- .cleanHub(data_S7[['DBhub']])

ColNames <- c('ADi','DBi','Confidence','ADORF','DBORF' )
data_S7[4:6] <- lapply(data_S7[4:6], function(x) {colnames(x) <- ColNames;x})
lapply(data_S7[4:6], head)
H2O2 <- data_S7[['H2O2']]       
MMS <- data_S7[['MMS']]
PoorCarbon <- data_S7[['PoorCarbon']]

Baseline <- data_S7[['Baseline']]
colnames(Baseline) <- Baseline[1,]
Baseline <- Baseline[-1, colnames(Baseline) %in% ColNames]

ListOfConditions <- list(Baseline, H2O2, MMS, PoorCarbon)
names(ListOfConditions) <- c('Baseline', 'H2O2', 'MMS', 'PoorCarbon')

##### add the PPi for each network
ListOfConditions <- lapply(ListOfConditions, function(x) {x$PPI <- .getPair(x); x})
lapply(ListOfConditions, head)

#saveRDS(ListOfConditions, 'Data/ListOfConditions_IS7.rds')
ListOfConditions <- readRDS('Data/ListOfConditions_IS7.rds')

## compare confidence distribution in different conditions
confDF <- rbind(data.frame(Conf=Baseline$Confidence, data='Baseline'),
                data.frame(Conf=H2O2$Confidence, data='H2O2'),
                data.frame(Conf=MMS$Confidence, data='MMS'),
                data.frame(Conf=PoorCarbon$Confidence, data='PoorCarbon'))
confDF$Conf <- as.numeric(confDF$Conf)
class(confDF$Conf)
p_conf <- ggplot(confDF ,aes(x=Conf, color= data))+geom_density()+
  theme_bw()+ggtitle('Confidence distribution')

ggplot(confDF ,aes(y=Conf,x=data , fill=data))+geom_boxplot()+
  theme_bw()+ggtitle('Confidence distribution')+scale_fill_brewer(palette="Set2")


