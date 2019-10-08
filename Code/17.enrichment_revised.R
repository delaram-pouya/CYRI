# https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/

### Function 
read.AND.clean <- function(File.path){
  data = read.csv(File.path, stringsAsFactors = F)
  data[,colnames(data) != 'X']}


add.global.Pvalues <- function(input.enriched, base_ORFS_localized, stress_ORFS_localized){
  test = merge(input.enriched, base_ORFS_localized, by.x='ORF', by.y='ORF', all.x=T)
  colnames(test)[ (ncol(input.enriched)+1) : ncol(test)] = paste0(colnames(test)[ (ncol(input.enriched)+1) : ncol(test)] , '_Baseline')
  
  test2 = merge(test, stress_ORFS_localized, by.x='ORF', by.y='ORF', all.x=T)
  colnames(test2)[ (ncol(test)+1) : ncol(test2)] = paste0(colnames(test2)[ (ncol(test)+1) : ncol(test2)] , '_stress')
  return(test2)
}


## import input data

listOfGlobal_GO_enrichments <- readRDS('Results/list_Of_GO_enriched_Global_localizations.rds')
base_ORFS_localized <- listOfGlobal_GO_enrichments$baseline
mms_ORFS_localized <- listOfGlobal_GO_enrichments$mms
h2o2_ORFS_localized <- listOfGlobal_GO_enrichments$h2o2
poorcarbon_ORFS_localized <- listOfGlobal_GO_enrichments$poorcarbon

stress_base_ORFS_list <- readRDS( 'Results/listOf_localPartenrs_stress_base.rds')
lapply(stress_base_ORFS_list, head)

mmsRes = read.AND.clean('Results/GO_enrich_tabs/mmsRes.csv')
h2o2Res = read.AND.clean('Results/GO_enrich_tabs/h2o2Res.csv')
poorcarbonRes = read.AND.clean('Results/GO_enrich_tabs/poorcarbonRes.csv')


##### adding the other 2 p values for that
mmsRes <- add.global.Pvalues(mmsRes, base_ORFS_localized, mms_ORFS_localized)
poorcarbonRes <- add.global.Pvalues(poorcarbonRes, base_ORFS_localized, poorcarbon_ORFS_localized)
h2o2Res <- add.global.Pvalues(h2o2Res, base_ORFS_localized, h2o2_ORFS_localized)

head(mmsRes)
head(poorcarbonRes)
head(h2o2Res)


## we are looking for cases in which the GO term is enriched between the baseline and stress condition,
## is also enriched in the stress condition but is not enriched in the baseline (or visa versa)
subset(test2, pValues <0.06)  ## YER034W
#### double check the code and make the slides for this section 
#### talk to DK abput the enrichemnt stuff
#### go over the candidates in the literature


##################################################################
## same result in the transpose 
input = matrix(c(12, 750, 86, 3253), 2, 2)
input.t = t(input)
fisher.test(input, alternative='two.sided')$p.value # greater
fisher.test(input.t, alternative='two.sided')$p.value # greater

factorial <- function(n) {
  if(n <= 1) {
    return(1)
  } else return(n * factorial(n-1))
}

n.choose.k <- function(n, k){
  factorial(n)/(factorial(k)*factorial(n-k))
}
n.choose.k(45, 11)

# Initialize variables
m <- 15       # Genes IN GO term
n <- 15       # Genes NOT IN GO term
k <- 15       # Gene hits, that is, differentially expressed
x <- c(0:15)  # Genes both IN GO term and differentially expressed 'hits'

# Use the dhyper built-in function for hypergeometric density
probabilities <- dhyper(x, m, n, k, log = FALSE)
probabilities
# Calculate the one-sided p-value for 12 or more genes both DE and IN GO term.
pvalue <- sum(probabilities[13:16])
pvalue
# Bar plot
data <- data.frame( x = x, y = probabilities )
ggplot(data, aes(x=factor(x), y=y)) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),
        axis.title.x=element_text(margin=margin(20,0,0,0)),
        axis.title.y=element_text(margin=margin(0,20,0,0))
  ) +
  geom_bar(stat="identity", fill=ifelse(data$x < 12,
                                        rgb(52, 73, 94, maxColorValue=255),
                                        rgb(231, 76, 60, maxColorValue=255)),
           colour="black") +
  labs(x = "DE genes IN GO term", y = "Probability")

