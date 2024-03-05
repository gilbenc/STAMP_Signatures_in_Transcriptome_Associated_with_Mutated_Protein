library("MetaGxBreast")

# get 39 optional studie names:
hub = ExperimentHub::ExperimentHub()
breastData <- query(hub, c("MetaGxBreast", "ExpressionSet"))
all_studies = breastData$title

# get METABRIC
esets = MetaGxBreast::loadBreastEsets(loadString = c("metabric"))
esets[[1]]

# get METABRIC expr. data
METABRIC_exp_data = esets$esets$METABRIC@assayData$exprs
METABRIC_exp_data = as.data.frame(METABRIC_exp_data)

rownames(METABRIC_exp_data[as.character(esets$esets$METABRIC@featureData@data$probeset),]) = as.character(esets$esets$METABRIC@featureData@data$gene)
PATH_data = "/home/gil/Desktop/GCN/datasets/Data_downloaded_directly_from_source/"
# write.csv(METABRIC_exp_data, paste0(PATH_data, "METABRIC_gene_expression_ILMN.csv"))

METABRIC_clinical_data = esets$esets$METABRIC@phenoData@data
# write.csv(METABRIC_clinical_data, paste0(PATH_data, "METABRIC_clinical_data.csv"))

esets = MetaGxBreast::loadBreastEsets(loadString = all_studies[!(all_studies %in% c("TCGA", "METABRIC"))])


GSES8644_exp_data = as.data.frame(esets$esets$GSE58644@assayData$exprs)
write.csv(GSES8644_exp_data, paste0(PATH_data, "GSES8644_exp_data.csv"))
# extract gene expression data from studies
all_studies_exp_data = as.data.frame(esets$esets$CAL@assayData$exprs)
for(study in names(esets$esets$GSE25066)){
  study_data = get(paste0("esets$esets$", study))
  exp_data = study_data@phenoData@data     
}


save.image("~/Desktop/t.RData")

esets_treatment = esets$esets[]

# Asses number of samples in each study
library(Biobase)
numSamples <- vapply(seq_along(esets$esets), FUN=function(i, esets){
 length(sampleNames(esets[[i]]))
 }, numeric(1), esets=esets$esets)
SampleNumberSummaryAll <- data.frame(NumberOfSamples = numSamples,
 row.names = names(esets$esets))
total <- sum(SampleNumberSummaryAll[,"NumberOfSamples"])
SampleNumberSummaryAll <- rbind(SampleNumberSummaryAll, total)
rownames(SampleNumberSummaryAll)[nrow(SampleNumberSummaryAll)] <- "Total"
require(xtable)
print(xtable(SampleNumberSummaryAll,digits = 2), floating = FALSE)

#pData Variables
pDataID <- c("er","pgr", "her2", "age_at_initial_pathologic_diagnosis",
                "grade", "dmfs_days", "dmfs_status", "days_to_tumor_recurrence",
                "recurrence_status", "days_to_death", "vital_status",
                "sample_type", "treatment")
 pDataPercentSummaryTable <- NULL
 pDataSummaryNumbersTable <- NULL
 pDataSummaryNumbersList <- lapply(esets$esets, function(x)
   vapply(pDataID, function(y) sum(!is.na(pData(x)[,y])), numeric(1)))

 pDataPercentSummaryList <- lapply(esets$esets, function(x)
   vapply(pDataID, function(y)
     sum(!is.na(pData(x)[,y]))/nrow(pData(x)), numeric(1))*100)
 pDataSummaryNumbersTable <- sapply(pDataSummaryNumbersList, function(x) x)
 pDataPercentSummaryTable <- sapply(pDataPercentSummaryList, function(x) x)
 rownames(pDataSummaryNumbersTable) <- pDataID
 rownames(pDataPercentSummaryTable) <- pDataID
 colnames(pDataSummaryNumbersTable) <- names(esets$esets)
 colnames(pDataPercentSummaryTable) <- names(esets$esets)
 pDataSummaryNumbersTable <- rbind(pDataSummaryNumbersTable, total)
 rownames(pDataSummaryNumbersTable)[nrow(pDataSummaryNumbersTable)] <- "Total"
 # Generate a heatmap representation of the pData
 pDataPercentSummaryTable <- t(pDataPercentSummaryTable)
 pDataPercentSummaryTable <- cbind(Name=(rownames(pDataPercentSummaryTable))
                                     ,pDataPercentSummaryTable)
 nba<-pDataPercentSummaryTable
 gradient_colors <- c("#ffffff","#ffffd9","#edf8b1","#c7e9b4","#7fcdbb",
                        "#41b6c4","#1d91c0","#225ea8","#253494","#081d58")
 library(lattice)
 nbamat<-as.matrix(nba)
 rownames(nbamat) <- nbamat[,1]
 nbamat <- nbamat[,-1]
 Interval <- as.numeric(c(10,20,30,40,50,60,70,80 ,90,100))
 levelplot(nbamat,col.regions=gradient_colors,
             main="Available Clinical Annotation",
             scales=list(x=list(rot=90, cex=1.5),
                           y= list(cex=1.5),key=list(cex=0.2)),
             at=seq(from=0,to=100,length=10),
             cex=0.2, ylab="", xlab="", lattice.options=list(),
             colorkey=list(at=as.numeric(factor(c(seq(from=0, to=100, by=10)))),
                             labels=as.character(c( "0","10%","20%","30%", "40%","50%",
                                                   "60%", "70%", "80%","90%", "100%"),
                                                   cex=0.2,font=1,col="brown",height=1,
                                                   width=1.4), col=(gradient_colors)))
 
 count = 0 
 for(i in 1:37) {
   count = count + sum(!is.na(esets$esets[[i]]$treatment))
 }
 print(count)
 