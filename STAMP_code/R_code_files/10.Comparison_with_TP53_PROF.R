### Analysis of GCN predictions through all cancer type models ###
library(ggplot2)
library(ggpubr)
library(dplyr)
library("RColorBrewer")

# see colors. change name to: Set1/Set2/Set3
# look here for more: http://www.sthda.com/english/wiki/colors-in-r#using-rcolorbrewer-palettes
display.brewer.pal(n = 8, name = 'Set3')
brewer.pal(n = 8, name = "Set3")
brw_set2 = brewer.pal(n = 8, name = "Set2")
data_type = "TP53_many_tumor_types"


PATH_GCN = paste0("//home//gil//Desktop//GCN//GCN_analysis//1.Prediction_tables//", data_type, "//GCN//")
PATH_ELR = paste0("//home//gil//Desktop//GCN//GCN_analysis//1.Prediction_tables//", data_type, "//ELR//")
if(data_type == "Pathways") PATH_plots <- "//home//gil//Desktop//GCN//GCN_analysis//Figures_for_paper//Mutations_TP53_PROF_comparison//TP53_Pathways_model//"
if(data_type == "TP53_many_tumor_types") PATH_plots <- "//home//gil//Desktop//GCN//GCN_analysis//Figures_for_paper//Mutations_TP53_PROF_comparison//"
analysis_PATH = "//home//gil//Desktop//GCN//GCN_analysis//"
datasets_PATH = "//home//gil//Desktop//GCN//datasets//Data_downloaded_directly_from_source//"

### Functions ###
# input: model, feature selection method.
# output: a dataframe of all samples with the respective predictions
create_table_for_model <- function(model, feature_selection=""){
  ### create one big pan cancer prediction table from all tumor specific GCN models.
  files <- list.files(paste0("//home//gil//Desktop//GCN//GCN_analysis//1.Prediction_tables//", data_type, "//", model, "//"))
  if(model == "GCN"){
    files <- files[grepl(files, pattern = "TP53_val")]  
  } else {
    files <- files[grepl(files, pattern = feature_selection)]
  }
  
  file <- files[1]
  predictions <- read.csv(paste0(get(paste0("PATH_", model)), file), row.names = "Sample.ID")
  predictions <- predictions[,-1]
  
  tumor_type_start_index <- gregexpr(pattern =paste0(model,"_"),file)[[1]] + 4
  tumor_type_end_index <- gregexpr(pattern ='_TP53',file)[[1]] - 1
  predictions$tumor_type <- substr(file, tumor_type_start_index, tumor_type_end_index)
  files <- files[-1]
  # read the rest of tumor types
  # create a file integrating GCN predictions from all tumor-specific models
  # and another file for the pan cancer predictions
  for(file in files) {
    predictions_one_tumor_type <- read.csv(paste0(get(paste0("PATH_", model)), file), row.names = "Sample.ID")
    predictions_one_tumor_type <- predictions_one_tumor_type[,-1]
    tumor_type_start_index <- gregexpr(pattern =paste0(model, '_'),file)[[1]] + 4
    tumor_type_end_index <- gregexpr(pattern ='_TP53',file)[[1]] - 1
    predictions_one_tumor_type$tumor_type <- substr(file, tumor_type_start_index, tumor_type_end_index)
    if(predictions_one_tumor_type$tumor_type[1] == "pan_cancer"){
      predictions_pan_cancer = predictions_one_tumor_type
    } else {
      predictions <- rbind(predictions, predictions_one_tumor_type)
    }
    
  }
  return(predictions)
  
}


# input: model, feature selection method.
# output: a dataframe of the pan cancer model
create_table_for_model_pan_cancer <- function(model, feature_selection=""){
  files <- list.files(paste0("//home//gil//Desktop//GCN//GCN_analysis//1.Prediction_tables//", data_type, "//", model, "//"))
  if(model == "GCN"){
    files <- files[grepl(files, pattern = "TP53_val")]  
  } else {
    files <- files[grepl(files, pattern = feature_selection)]
  }
  
  file <- files[grepl(files, pattern = "pan")]
  predictions <- read.csv(paste0(get(paste0("PATH_", model)), file), row.names = "Sample.ID")
  predictions <- predictions[,-1]
  
  tumor_type_start_index <- gregexpr(pattern =paste0(model,"_"),file)[[1]] + 4
  tumor_type_end_index <- gregexpr(pattern ='_TP53',file)[[1]] - 1
  predictions$tumor_type <- substr(file, tumor_type_start_index, tumor_type_end_index)
  return(predictions)
}


### Create Linear score
# create a function to calculate the logit of a value
logit <- function(x) {
  logit_value <- log(x/(1-x))
  if(logit_value==Inf){
    return(17)
  }
  if(logit_value == -Inf){
    return(-105)
  }
  return(logit_value)
}


expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  
  y <- unique(y)
  
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  
  do.call(rbind, lapply(seq_along(x), g))
}


### DATA ### 

# LOAD MODELS
# ELR_predictions <- create_table_for_model("ELR", "mutsigdb")
GCN_predictions <- create_table_for_model("GCN")

# LOAD cbioportal mutations data
# get protein variants for TCGA sample IDs
TCGA_TP53_mutated_pan_cancer <- read.csv(paste0(datasets_PATH, "cbioportal_TCGA_pan_cancer_TP53_mutated_samples1.0.tsv"), sep = "\t")
TCGA_TP53_mutated_pan_cancer <- TCGA_TP53_mutated_pan_cancer[,c("Sample.ID", "Protein.Change", "Mutation.Type", "Annotation", "Functional.Impact")]

# LOAD TP53_PROF predictions
TP53_PROF_preds <- read.csv(paste0(datasets_PATH, "TP53_PROF_variant_predictions.csv"))
TP53_PROF_preds <- TP53_PROF_preds[,c("protein_change", "label","predictions", "score")]
colnames(TP53_PROF_preds) <- c("Protein.Change", "TP53_PROF_label", "TP53_PROF_predictions", "TP53_PROF_score")
TP53_PROF_preds$Protein.Change <- sapply(TP53_PROF_preds$Protein.Change, function(x){substr(x, 3, 7)})

# LOAD TP53_PROF with features
TP53_PROF_preds_with_features <- read.csv(paste0(datasets_PATH, "TP53_PROF_variant_predictions_with_features.csv"))
TP53_PROF_preds_with_features <- TP53_PROF_preds_with_features[,c("Giac_A549_WT_Nut_norm", "Kolt_RFS_H1299_norm", "WAF1_percent", "percent_14_3_3_s", "protein_change", "label","predictions", "score")]
colnames(TP53_PROF_preds_with_features) <- c("Giac_A549_WT_Nut_norm", "Kolt_RFS_H1299_norm", "WAF1_percent", "percent_14_3_3_s","Protein.Change", "TP53_PROF_label", "TP53_PROF_predictions", "TP53_PROF_score")
TP53_PROF_preds_with_features$Protein.Change <- sapply(TP53_PROF_preds_with_features$Protein.Change, function(x){substr(x, 3, 7)})




### CODE ###

GCN_predictions$linear_score <- sapply(GCN_predictions$pred_score, logit)
GCN_predictions$Sample.ID <- row.names(GCN_predictions)


# save table with predictions and linear score
# write.csv(GCN_predictions, paste0(analysis_PATH, "1.Prediction_tables//GCN_predictions_linear_score_tumor_specific_approach_24_cancer_types.csv"))

### 4. Create discrepancies table
GCN_discrepancies <- GCN_predictions[GCN_predictions$label != GCN_predictions$prediction,]

### 5. Test discrepancies with model's predictions

# unique(TCGA_TP53_mutated_pan_cancer$Mutation.Type)

GCN_preds_with_protein_variant_all_mutations <- merge(GCN_predictions, TCGA_TP53_mutated_pan_cancer, "Sample.ID")
GCN_preds_with_protein_variant = GCN_preds_with_protein_variant_all_mutations[GCN_preds_with_protein_variant_all_mutations$Mutation.Type=="Missense_Mutation",]

# Modify mutation types for this analysis (similar to Thierry's analysis in https://www.nature.com/articles/s41598-020-74892-2 )
GCN_preds_with_protein_variant_all_mutations$Mutation.Type[GCN_preds_with_protein_variant_all_mutations$Mutation.Type %in% 
                                                             c("Frame_Shift_Del", "Frame_Shift_Ins")] <- "Frame_Shift"
GCN_preds_with_protein_variant_all_mutations$Mutation.Type[GCN_preds_with_protein_variant_all_mutations$Mutation.Type %in% 
                                                             c("Fusion", "Nonsense_Mutation")] <- "Nonsense_Mutation"
GCN_preds_with_protein_variant_all_mutations$Mutation.Type[GCN_preds_with_protein_variant_all_mutations$Mutation.Type %in% 
                                                             c("In_Frame_Del", "In_Frame_Ins")] <- "In_Frame"
GCN_preds_with_protein_variant_all_mutations$Mutation.Type[GCN_preds_with_protein_variant_all_mutations$Mutation.Type %in% 
                                                             c("Splice_Region", "Splice_Site", "Translation_Start_Site")] <- "Splice_Site"

table(GCN_preds_with_protein_variant_all_mutations$Mutation.Type)

# comparisons used in boxplot- all mutation types
all_comparisons = expand.grid.unique(GCN_preds_with_protein_variant_all_mutations$Mutation.Type,GCN_preds_with_protein_variant_all_mutations$Mutation.Type)
all_comparisons = all_comparisons[all_comparisons[,1]!= "Translation_Start_Site" & all_comparisons[,2]!= "Translation_Start_Site",]

my_comparisons = mapply(c, all_comparisons[,1], all_comparisons[,2], SIMPLIFY = FALSE)



### TP53_PROF analysis of missense variants in different tumor types.

GCN_preds_with_variant <- merge(GCN_preds_with_protein_variant, TP53_PROF_preds, "Protein.Change")

GCN_preds_with_protein_variant <- merge(GCN_preds_with_protein_variant, TP53_PROF_preds_with_features, "Protein.Change")



### PAN CANCER MODEL ###
# Run same code for pan cancer model
GCN_predictions_pan_cancer <- create_table_for_model_pan_cancer("GCN")
GCN_predictions_pan_cancer$linear_score <- sapply(GCN_predictions_pan_cancer$pred_score, logit)
GCN_predictions_pan_cancer$Sample.ID <- row.names(GCN_predictions_pan_cancer)
GCN_discrepancies_pan_cancer <- GCN_predictions_pan_cancer[GCN_predictions_pan_cancer$label != GCN_predictions_pan_cancer$prediction,]
GCN_preds_with_protein_variant_all_mutations_pan_cancer <- merge(GCN_predictions_pan_cancer, TCGA_TP53_mutated_pan_cancer, "Sample.ID")
GCN_preds_with_protein_variant_pan_cancer = GCN_preds_with_protein_variant_all_mutations_pan_cancer[GCN_preds_with_protein_variant_all_mutations_pan_cancer$Mutation.Type=="Missense_Mutation",]
# Modify mutation types for this analysis (similar to Thierry's analysis in https://www.nature.com/articles/s41598-020-74892-2 )
GCN_preds_with_protein_variant_all_mutations_pan_cancer$Mutation.Type[GCN_preds_with_protein_variant_all_mutations_pan_cancer$Mutation.Type %in% 
                                                                        c("Frame_Shift_Del", "Frame_Shift_Ins")] <- "Frame_Shift"
GCN_preds_with_protein_variant_all_mutations_pan_cancer$Mutation.Type[GCN_preds_with_protein_variant_all_mutations_pan_cancer$Mutation.Type %in% 
                                                                        c("Fusion", "Nonsense_Mutation")] <- "Nonsense_Mutation"
GCN_preds_with_protein_variant_all_mutations_pan_cancer$Mutation.Type[GCN_preds_with_protein_variant_all_mutations_pan_cancer$Mutation.Type %in% 
                                                                        c("In_Frame_Del", "In_Frame_Ins")] <- "In_Frame"
GCN_preds_with_protein_variant_all_mutations_pan_cancer$Mutation.Type[GCN_preds_with_protein_variant_all_mutations_pan_cancer$Mutation.Type %in% 
                                                                        c("Splice_Region", "Splice_Site", "Translation_Start_Site")] <- "Splice_Site"
table(GCN_preds_with_protein_variant_all_mutations_pan_cancer$Mutation.Type)
all_comparisons = expand.grid.unique(GCN_preds_with_protein_variant_all_mutations_pan_cancer$Mutation.Type,GCN_preds_with_protein_variant_all_mutations_pan_cancer$Mutation.Type)
all_comparisons = all_comparisons[all_comparisons[,1]!= "Translation_Start_Site" & all_comparisons[,2]!= "Translation_Start_Site",]
my_comparisons = mapply(c, all_comparisons[,1], all_comparisons[,2], SIMPLIFY = FALSE)
GCN_preds_with_variant_pan_cancer <- merge(GCN_preds_with_protein_variant_pan_cancer, TP53_PROF_preds, "Protein.Change")
GCN_preds_with_protein_variant_pan_cancer <- merge(GCN_preds_with_protein_variant_pan_cancer, TP53_PROF_preds_with_features, "Protein.Change")
GCN_preds_with_protein_variant_pan_cancer$TP53_PROF_label[GCN_preds_with_protein_variant_pan_cancer$TP53_PROF_label == "null"] <- 
                                                                        GCN_preds_with_protein_variant_pan_cancer$TP53_PROF_predictions[
                                                                          GCN_preds_with_protein_variant_pan_cancer$TP53_PROF_label == "null"]
GCN_discrepancies_with_protein_variant_pan_cancer <- GCN_preds_with_protein_variant_pan_cancer[
  GCN_preds_with_protein_variant_pan_cancer$label != GCN_preds_with_protein_variant_pan_cancer$prediction,]
# See intersection of GCN preds and TP53_PROF preds in *** Missense mutated samples ***
table(GCN_preds_with_protein_variant_pan_cancer$prediction, GCN_preds_with_protein_variant_pan_cancer$TP53_PROF_label)
# chi sq. between PROF and GNC preds:
chisq.test(GCN_preds_with_protein_variant_pan_cancer$prediction, GCN_preds_with_protein_variant_pan_cancer$TP53_PROF_label)
# t.test for linear score and PROF
t.test(GCN_preds_with_protein_variant_pan_cancer$TP53_PROF_label=="D", GCN_preds_with_protein_variant_pan_cancer$linear_score)
wilc_pan_cancer = wilcox.test(GCN_preds_with_protein_variant_pan_cancer$linear_score~GCN_preds_with_protein_variant_pan_cancer$TP53_PROF_label,
                              alternative = "greater")
# print(wilc_tumor_type)
tiff(paste0(PATH_plots, "GCN_TP53_PROF_comparison_TP53mutated_samples_boxplot_PAN_CANCER_MODEL_wilc_p_",
            round(wilc_pan_cancer$p.value,5),".tiff"), units="in", width=5, height=4, res=300)
print(ggplot(GCN_preds_with_protein_variant_pan_cancer, aes(x=TP53_PROF_label, y=linear_score, fill=TP53_PROF_label)) +
        geom_boxplot()+ labs(x = paste0("TP53_PROF prediction\n Pan Cancer"), y = "GCN Linear score", fill = "TP53_PROF") + scale_fill_manual(values = c("#ef8a62", "#67a9cf")) +
        stat_compare_means(data = GCN_preds_with_protein_variant_pan_cancer, method = "wilcox.test",label.x=1, label.y=min(GCN_preds_with_protein_variant_pan_cancer$linear_score),
                           method.args = list(alternative = "less"), paired = F, size = 6) + theme_bw()+geom_jitter(color=brw_set2[3], size=0.4, alpha=0.8)+
        font("x.text", size = 20) +font("y.text", size = 20) + font("xy.title", size = 20) + 
        theme(strip.text.x = element_text(size=20))+ theme(legend.text=element_text(size=20))#+ theme(legend.title=element_text(size=15))
        )
dev.off()
########## FOR RESULTS: notice table for pan cancer mutated samples. 74 ND variants. ***35*** of them are in the discrepancies, i.e. predicted ND by TP53_PROF despite having a mutation.
table(GCN_preds_with_protein_variant_pan_cancer$label, GCN_preds_with_protein_variant_pan_cancer$TP53_PROF_label)
table(GCN_discrepancies_with_protein_variant_pan_cancer$label, GCN_discrepancies_with_protein_variant_pan_cancer$TP53_PROF_label)
### Mutation types analysis
my_comparisons = list(c("Frame Shift", "Missense ND"), c("Nonsense", "Missense ND"), c("Splice", "Missense ND"),c("Missense D", "Missense ND"), c("In Frame", "Missense ND"))
for(mut in unique(GCN_preds_with_protein_variant_pan_cancer$Protein.Change)){
  if("ND" %in% GCN_preds_with_protein_variant_pan_cancer$TP53_PROF_label[GCN_preds_with_protein_variant_pan_cancer$Protein.Change==mut]){
    GCN_preds_with_protein_variant_all_mutations_pan_cancer$Mutation.Type[GCN_preds_with_protein_variant_all_mutations_pan_cancer$Protein.Change == mut] = "ND"
  } else {
    GCN_preds_with_protein_variant_all_mutations_pan_cancer$Mutation.Type[GCN_preds_with_protein_variant_all_mutations_pan_cancer$Protein.Change == mut] = "D"
  }
}
# GCN_preds_with_protein_variant_all_mutations_pan_cancer = GCN_preds_with_protein_variant_all_mutations_pan_cancer[GCN_preds_with_protein_variant_all_mutations_pan_cancer$Mutation.Type!= "In_Frame",]
GCN_preds_with_protein_variant_all_mutations_pan_cancer$Mutation.Type = factor(GCN_preds_with_protein_variant_all_mutations_pan_cancer$Mutation.Type, 
                                                                    levels=c("Frame_Shift","Nonsense_Mutation", "Splice_Site", "D", "ND", "In_Frame"))
levels(GCN_preds_with_protein_variant_all_mutations_pan_cancer$Mutation.Type) = c("Frame Shift", "Nonsense", "Splice", "Missense D", "Missense ND", "In Frame")
tiff(paste0(PATH_plots, "boxplot_GCN_linear_Mutation_type_PAN_CANCER_MODEL.tiff"), units="in", width=18, height=5.5, res=300)
ggboxplot(GCN_preds_with_protein_variant_all_mutations_pan_cancer, x= "Mutation.Type", y = "linear_score",
          color = "Mutation.Type")+labs(x = paste0("Mutation Type\nPan Cancer"), y = "GCN Linear score", fill = "Mutation Type")+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "greater"), size=6)+
  geom_jitter(color=brw_set2[3], size=0.4, alpha=0.8)+font("x.text", size = 20) +font("y.text", size = 20) + font("xy.title", size = 20) + 
  theme(strip.text.x = element_text(size=20))+theme(legend.position="none")#+ theme(legend.text=element_text(size=20))#+ theme(legend.title=element_text(size=20))
dev.off()


# boxplot mutations pan cancer
my_comparisons = list(c("Frame Shift", "Missense ND"), c("Nonsense", "Missense ND"), c("Splice", "Missense ND"),c("Missense D", "Missense ND"), c("Missense D", "In Frame"), c("Frame Shift", "In Frame"), c("Nonsense", "In Frame"), c("Splice", "In Frame"))




GCN_preds_with_protein_variant$TP53_PROF_label[GCN_preds_with_protein_variant$TP53_PROF_label == 
                                                 "null"] <- GCN_preds_with_protein_variant$TP53_PROF_predictions[
                                                   GCN_preds_with_protein_variant$TP53_PROF_label == "null"]


GCN_discrepancies_with_protein_variant <- GCN_preds_with_protein_variant[
  GCN_preds_with_protein_variant$label != GCN_preds_with_protein_variant$prediction,]




# See intersection of GCN preds and TP53_PROF preds in *** Missense mutated samples ***
table(GCN_preds_with_protein_variant$prediction, GCN_preds_with_protein_variant$TP53_PROF_label)
# chi sq. between PROF and GNC preds:
chisq.test(GCN_preds_with_protein_variant$prediction, GCN_preds_with_protein_variant$TP53_PROF_label)
# t.test for linear score and PROF
t.test(GCN_preds_with_protein_variant$TP53_PROF_label=="D", GCN_preds_with_protein_variant$linear_score)

display.brewer.pal(n = 5, name = 'PuBu')
brw_dark2 <- brewer.pal(n = 3, name = "Dark2")
brw_set2 <- brewer.pal(n = 5, name = "Set2")
library(ggpubr)

# Examine PROF predictions in specific tumor types
for(tumor_type in unique(GCN_preds_with_protein_variant$tumor_type)){
  GCN_preds_with_protein_variant_tumor_type = GCN_preds_with_protein_variant[GCN_preds_with_protein_variant$tumor_type==tumor_type,]
  
  if(nrow(GCN_preds_with_protein_variant_tumor_type)>50 & 
       sum(GCN_preds_with_protein_variant_tumor_type$TP53_PROF_label=="ND") > 3){
      # Show tumor type, table and fit wilcoxon test
      # print(tumor_type)
      # print(table(GCN_preds_with_protein_variant_tumor_type$TP53_PROF_label))
      wilc_tumor_type = wilcox.test(GCN_preds_with_protein_variant_tumor_type$linear_score~GCN_preds_with_protein_variant_tumor_type$TP53_PROF_label,
                                    alternative = "greater")
      # print(wilc_tumor_type)
      

      
      tiff(paste0(PATH_plots, "GCN_TP53_PROF_comparison_TP53mutated_samples_boxplot_",tumor_type,"_wilc_p_",
                  round(wilc_tumor_type$p.value,5),".tiff"), units="in", width=5, height=4, res=300)
      print(ggplot(GCN_preds_with_protein_variant_tumor_type, aes(x=TP53_PROF_label, y=linear_score, fill=TP53_PROF_label)) +
        geom_boxplot()+ labs(x = paste0("TP53_PROF prediction\n",tumor_type), y = "GCN Linear score", fill = "TP53_PROF") + scale_fill_manual(values = c("#ef8a62", "#67a9cf")) +
        stat_compare_means(data = GCN_preds_with_protein_variant_tumor_type, method = "wilcox.test",label.x=1.2, label.y=min(GCN_preds_with_protein_variant_tumor_type$linear_score),
                          method.args = list(alternative = "less"), paired = F, size = 6) +geom_jitter(color=brw_set2[3], size=0.4, alpha=0.8)+
        theme_bw()+font("x.text", size = 20) +font("y.text", size = 20) + font("xy.title", size = 20) + 
          theme(strip.text.x = element_text(size=20))+ theme(legend.text=element_text(size=20)))#+ theme(legend.title=element_text(size=20)))
      dev.off()
  }
}



# remove outliers
GCN_preds_with_protein_variant_no_outliers = GCN_preds_with_protein_variant[GCN_preds_with_protein_variant$linear_score > -31,]

wilc_pan_cancer = wilcox.test(GCN_preds_with_protein_variant_no_outliers$linear_score~GCN_preds_with_protein_variant_no_outliers$TP53_PROF_label,
                              alternative = "greater")
print(wilc_pan_cancer)

PATH_plots <- "//home//gil//Desktop//GCN//GCN_analysis//Figures_for_paper//Mutations_TP53_PROF_comparison//"

tiff(paste0(PATH_plots, "GCN_TP53_PROF_comparison_TP53mutated_samples_boxplot_pan_cancer_wilc_p_",round(wilc_pan_cancer$p.value,4),".tiff"), units="in", width=5, height=3, res=300)
ggplot(GCN_preds_with_protein_variant_no_outliers, aes(x=TP53_PROF_label, y=linear_score, fill=TP53_PROF_label)) +
  geom_boxplot()+ labs(x = "TP53_PROF prediction\nPan Cancer", y = "GCN Linear score", fill = "TP53_PROF") + scale_fill_manual(values = c(brw_set2[5], brw_set2[4])) +
  stat_compare_means(data = GCN_preds_with_protein_variant_no_outliers, method = "wilcox.test", label.x=1.3, label.y=-30,
                     method.args = list(alternative = "less"), paired = F, size = 3.5) +geom_jitter(color=brw_set2[3], size=0.4, alpha=0.8)+
  theme_bw()
dev.off()

library("gplots")

GCN_preds_with_protein_variant$GCN_norm_pred = as.numeric(
                                            GCN_preds_with_protein_variant$linear_score>mean(GCN_preds_with_protein_variant$linear_score)
                                            )
# Chi sqaure and balloon plots
table_for_PROF_GCN = table(GCN_preds_with_protein_variant$prediction, GCN_preds_with_protein_variant$TP53_PROF_label)
cur_chisq = chisq.test(table_for_PROF_GCN)
print(cur_chisq)
cur_chisq$expected

# the balooneplot is an intuitive visualisation tool that quantifies integers as "baloons"
balloonplot(table_for_PROF_GCN, main ="GCN and TP53_PROF preds", xlab ="", ylab="",label = FALSE, show.margins = FALSE)
as.table(cur_chisq$expected)

# now we show a ballonplot with the difference between our real data and the expected data
balloonplot(table_for_PROF_GCN - as.table(cur_chisq$expected), main ="GCN and TP53_PROF pred", xlab ="", ylab="",label = FALSE, show.margins = FALSE)


### Generate the pan cancer analysis

for(mut in unique(GCN_preds_with_protein_variant$Protein.Change)){
  if("ND" %in% GCN_preds_with_protein_variant$TP53_PROF_label[GCN_preds_with_protein_variant$Protein.Change==mut]){
    GCN_preds_with_protein_variant_all_mutations$Mutation.Type[GCN_preds_with_protein_variant_all_mutations$Protein.Change == mut] = "ND"
  } else {
    GCN_preds_with_protein_variant_all_mutations$Mutation.Type[GCN_preds_with_protein_variant_all_mutations$Protein.Change == mut] = "D"
  }
}


GCN_preds_with_protein_variant_all_mutations$Mutation.Type = factor(GCN_preds_with_protein_variant_all_mutations$Mutation.Type, 
                                                                    levels=c("Frame_Shift","Nonsense_Mutation", "Splice_Site", "D", "ND", "In_Frame" ))
levels(GCN_preds_with_protein_variant_all_mutations$Mutation.Type) = c("Frame Shift", "Nonsense", "Splice", "Missense D", "Missense ND", "In Frame")

# boxplot mutations pan cancer
my_comparisons = list(c("Frame Shift", "Missense ND"), c("Nonsense", "Missense ND"), c("Splice", "Missense ND"),c("Missense D", "Missense ND"), c("Missense D", "In Frame"))#, c("Frame Shift", "In Frame"), c("Nonsense", "In Frame"), c("Splice", "In Frame"))



tiff(paste0(PATH_plots, "boxplot_GCN_linear_Mutation_type_pan_cancer.tiff"), units="in", width=9, height=5.5, res=300)
ggboxplot(GCN_preds_with_protein_variant_all_mutations, x= "Mutation.Type", y = "linear_score",
          color = "Mutation.Type")+labs(x = paste0("Mutation Type\nPan Cancer"), y = "GCN Linear score", fill = "Mutation Type")+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "greater"), size=4.5)
dev.off()

for(tumor_type in c("BRCA", "LUAD")){#unique(GCN_preds_with_protein_variant$tumor_type)){
  if(nrow(GCN_preds_with_protein_variant_all_mutations[GCN_preds_with_protein_variant_all_mutations$tumor_type==tumor_type,])>50 &
     sum(GCN_preds_with_protein_variant_all_mutations[GCN_preds_with_protein_variant_all_mutations$tumor_type==tumor_type,"Mutation.Type"]=="Missense ND") > 3) {
    
    tiff(paste0(PATH_plots, "boxplot_GCN_linear_Mutation_type_",tumor_type,"_", data_type,".tiff"), units="in", width=12, height=7, res=300)
    print(ggboxplot(GCN_preds_with_protein_variant_all_mutations[GCN_preds_with_protein_variant_all_mutations$tumor_type==tumor_type,], x= "Mutation.Type", y = "linear_score",
                    color = "Mutation.Type")+labs(x = paste0("Mutation Type\n",tumor_type), y = "GCN Linear score", fill = "Mutation Type")+
            stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "greater"), size=7)+
            font("x.text", size = 20) +font("y.text", size = 20) + font("xy.title", size = 20) + #geom_jitter(color=brw_set2[3], size=0.4, alpha=0.8)+
            theme(strip.text.x = element_text(size=20))+ theme(legend.text=element_text(size=20))+ theme(legend.title=element_text(size=20)))
    dev.off()
  }
}


