library(org.Sc.sgd.db)

###### generate list of all GO annotations
all_GO_terms <- unique(c(unlist(lapply(ListOfConditions, function(x) x$ADORF)), 
                         unlist(lapply(ListOfConditions, function(x) x$DBORF))))
all_GO_terms <- data.frame(all_GO_terms)
head(all_GO_terms)
write.table(all_GO_terms, 'Data/GO/all_GO_terms.txt', quote = F, row.names = F, col.names = F)

##### upload them online to get the Gene annotation results , and import them 
GOterms <- read.delim('Data/GO/GO_terms_results.tsv', header = F)
colnames(GOterms) <- c('SGD_id', 'sys_name', 'Organism', 'stand_name', 'Description')
head(GOterms)
sum(is.na(GOterms))

###### GO slim mapping to map genes to pathways
GOmap <- read.delim('Data/GO/go_slim_mapping.tab', header = F)
ListOfConditions <- readRDS('Data/ListOfConditions_IS7.rds')
colnames(GOmap)<- c('sys_name','stand_name','SGD_id', 'CFP','Pathway', 'GOnumber', 'type')

GOdata <- merge(GOterms, GOmap, by.x='SGD_id', by.y='SGD_id', all.x= T)
GOdata[is.na(GOdata)] <- ''

sum(GOdata$sys_name.x != GOdata$sys_name.y)
GOdata[GOdata$sys_name.x != GOdata$sys_name.y, ]

head(GOdata)
saveRDS(GOdata, 'Data/GO/GO_terms_pathway.rds')
GOdata <- readRDS('Data/GO/GO_terms_pathway.rds')




test <- read.delim('Data/subcellular_location.tsv')
head(test)
table(test$Main.location)
table(unlist(str_split(test$Main.location,';')))

