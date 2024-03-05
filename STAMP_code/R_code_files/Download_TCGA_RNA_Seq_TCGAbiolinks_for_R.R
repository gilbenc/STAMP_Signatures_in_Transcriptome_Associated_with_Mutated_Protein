###########################################
### Code to download TCGa rna-seq data  ###
### and mutations per gene              ###
###########################################

# use this library to install TCGAbiolinks
# library(BiocManager)
library(TCGAbiolinks)



############
### Main ###
############


### Part 1: expression
exp.path <- "E://GDCdata//Expression_tables/"

 
file_type = "normalized_results"
to.download = F
projects <- TCGAbiolinks::getGDCprojects()$project_id
projects <- projects[grepl("TCGA", projects)]

expression_data <- NULL
sum <- 0
for(i in 1:length(projects)) {
        query <- GDCquery(project = projects[i],
                          data.category = "Gene expression",
                          data.type = "Gene expression quantification",
                          platform = "Illumina HiSeq",
                          file.type = file_type,
                          experimental.strategy = "RNA-Seq",
                          legacy = TRUE)
        
        if(to.download) {
                GDCdownload(query, directory = "E://GDCdata//", method = "api", files.per.chunk =10)
        }
        # # data will contain expression data for TCGA-BRCA project
        data_norm <- GDCprepare(query,
                                save=TRUE,
                                directory = "E://GDCdata/",
                                save.filename = "Gene_Expression_Quantification.rda",
                                summarizedExperiment = FALSE)
        
        data_norm <- as.data.frame(t(data_norm))
        data_norm$project <- projects[i]
        expression_data <- rbind(expression_data, data_norm)
        sum = sum + nrow(data_norm)
        print(paste0("finished with ", projects[i], ", with ", nrow(data_norm), " samples. total samples so far: ", sum, ", project number: ", i, " of ", length(projects), " projects."))
}

write.csv(expression_data, "E://GDCdata//Expression_tables//RNASeq_33_TCGA_projects_downloaded_with_TCGAbiolinks_2021_03_01.csv")


### Part 2: mutations ###
mut.path <- "E://GDCdata//Mutations_tables_from_TCGAbiolinks//"


# get all profiled samples per cancer_type in TCGA and label them by gene's mutational status.
# input: gene (hugo symbol) and should a csv be written to disc?
get_mutation_profile_from_TCGAbiolinks <- function(gene, csv_save) {
        mutations_data <- NULL
        
        for(i in 1:length(projects)) {
                tumor = strsplit(projects[i], "-")[[1]][2]
                query <- TCGAbiolinks::GDCquery_Maf(tumor = tumor, pipelines = "mutect2")
                # get all sample IDs that were queried for this cancer type
                Sample_IDs_profile <- substr(unique(query$Tumor_Sample_Barcode), 1, 15)
                df = data.frame(Sample_IDs_profile, tumor, 0)
                colnames(df) <- c("Sample ID", "tumor_type", gene)
                # get only query values relevant for given gene
                query <- query[query$Hugo_Symbol == gene,]
                Sample_IDs_mutated <- substr(unique(query$Tumor_Sample_Barcode), 1, 15)
                if(length(Sample_IDs_mutated)>0)
                        df[df$`Sample ID` %in% Sample_IDs_mutated, gene] <- 1
                mutations_data <- rbind(mutations_data, df)        
                
                print(paste0(i, ": ", projects[i], " had ", length(Sample_IDs_profile), " samples profiled and ", length(Sample_IDs_mutated), " mutated in ", gene))
        }
        
        mutations_data <- read.csv(paste0(mut.path, gene, "_Mutations_from_TCGAbiolinks.csv"))
        
        mutations_data <- mutations_data[!duplicated(mutations_data),]
        projects_mut_df <- data.frame(unique(mutations_data$tumor) ,sapply(unique(mutations_data$tumor), function(x){sum(mutations_data$tumor == x)}), sapply(unique(mutations_data$tumor), function(x){sum(mutations_data[mutations_data$tumor == x, gene])}))
        colnames(projects_mut_df) <- c("project", "profiled samples", paste0(gene, " mutated samples"))
        
        if(csv_save) {
                write.csv(x = mutations_data, row.names = F, file = paste0(mut.path, gene, "_Mutations_from_TCGAbiolinks.csv"))
                write.csv(x = projects_mut_df, row.names = F, file = paste0(mut.path, gene, "_Mutations_from_TCGAbiolinks_documentation.csv"))    
        }
        
}

# gene's hugo symbol
gene = "TP53"
csv_save = T
get_mutation_profile_from_TCGAbiolinks(gene, csv_save)


### Combine TCGAbiolinks and cbioportal profiled sample IDs ###
TCGAbiolinks_mutations_data <- mutations_data
# load cbioportal datasets
cbioportal_clinical_data <- read.csv("E://GDCdata//cBioPortal_tables_Mutations_clinical//cbioportal_TCGApancanceratlas_mutations_profiled_samples_clinical_data.tsv", sep = "\t")
cbioportal_clinical_data <- cbioportal_clinical_data[c("Sample.ID", "TCGA.PanCanAtlas.Cancer.Type.Acronym")]
cbioportal_mutations_data <- read.csv("E://GDCdata//cBioPortal_tables_Mutations_clinical//cbioportal_TCGApancanceratlas_mutations_profiled_samples_TP53_mutated.tsv", sep = "\t")
cbioportal_mutations_data <- cbioportal_mutations_data[,c("Sample.ID", "Cancer.Type")]

# create label for cbioportal samples
cbioportal_clinical_data$TP53 = 0
cbioportal_clinical_data$TP53[cbioportal_clinical_data$Sample.ID %in% cbioportal_mutations_data$Sample.ID] <- 1
colnames(cbioportal_clinical_data)[2] <- "tumor_type"

# create combined dataset with labels
TCGAbiolinks_cbioportal_combined_profiled_sampleIDs <- TCGAbiolinks_mutations_data
cbioportal_clinical_data_not_in_TCGAbiolinks <- cbioportal_clinical_data[!cbioportal_clinical_data$Sample.ID %in% TCGAbiolinks_mutations_data$Sample.ID,]
TCGAbiolinks_cbioportal_combined_profiled_sampleIDs <- rbind(TCGAbiolinks_mutations_data, cbioportal_clinical_data_not_in_TCGAbiolinks)

# 
projects_mut_df <- data.frame(unique(TCGAbiolinks_cbioportal_combined_profiled_sampleIDs$tumor_type) ,sapply(unique(TCGAbiolinks_cbioportal_combined_profiled_sampleIDs$tumor_type), function(x){sum(TCGAbiolinks_cbioportal_combined_profiled_sampleIDs$tumor_type == x)}), sapply(unique(TCGAbiolinks_cbioportal_combined_profiled_sampleIDs$tumor_type), function(x){sum(TCGAbiolinks_cbioportal_combined_profiled_sampleIDs[TCGAbiolinks_cbioportal_combined_profiled_sampleIDs$tumor_type == x, gene])}))
colnames(projects_mut_df) <- c("project", "profiled samples", paste0(gene, " mutated samples"))

# write combined data down
write.csv(x = TCGAbiolinks_cbioportal_combined_profiled_sampleIDs, file ="E://GDCdata//TCGAbiolinks_cbioportal_combined_profiled_sample_IDs_with_labels.csv", row.names = F)
write.csv(x = projects_mut_df, row.names = F, file= "E://GDCdata//TCGAbiolinks_cbioportal_combined_projects_summary.csv")

