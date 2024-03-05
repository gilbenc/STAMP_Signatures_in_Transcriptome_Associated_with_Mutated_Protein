### ROC curves for the different models ###
library(ggplot2)
library(ggpubr)
library("RColorBrewer")
library(ROCR)
library(tidyr)
library(pROC)

### Functions ###
create_table_for_model <- function(data_type, model, feature_selection, pan){
  assign(paste0("PATH_", model), 
         paste0("//home//gil//Desktop//GCN//GCN_analysis//1.Prediction_tables//", data_type, "//", model,"//"))
  ### 1. create one big pan cancer prediction table from all tumor specific GCN models.
  files <- list.files(paste0("//home//gil//Desktop//GCN//GCN_analysis//1.Prediction_tables//", data_type, "//", model, "//"))
  if(model == "GCN"){
    files <- files[grepl(files, pattern = "TP53_val_AUC")]  
  } else {
    files <- files[grepl(files, pattern = feature_selection)]
  }
  
  file <- files[1]
  predictions <- read.csv(paste0(get(paste0("PATH_", model)), file), row.names = "Sample.ID")
  predictions <- predictions[,-1]
  
  tumor_type_start_index <- gregexpr(pattern =paste0(model, '_'),file)[[1]] + 1 + nchar(model)
  tumor_type_end_index <- gregexpr(pattern ='_TP53',file)[[1]] - 1
  predictions$tumor_type <- substr(file, tumor_type_start_index, tumor_type_end_index)
  files <- files[-1]
  # read the rest of tumor types
  # create a file integrating GCN predictions from all tumor-specific models
  # and another file for the pan cancer predictions
  for(file in files) {
    predictions_one_tumor_type <- read.csv(paste0(get(paste0("PATH_", model)), file), row.names = "Sample.ID")
    predictions_one_tumor_type <- predictions_one_tumor_type[,-1]
    tumor_type_start_index <- gregexpr(pattern =paste0(model, '_'),file)[[1]] + 1 + nchar(model)
    tumor_type_end_index <- gregexpr(pattern ='_TP53',file)[[1]] - 1
    predictions_one_tumor_type$tumor_type <- substr(file, tumor_type_start_index, tumor_type_end_index)
    # if(predictions_one_tumor_type$tumor_type[1] == "pan_cancer"){
    #   predictions_pan_cancer = predictions_one_tumor_type
    # } else {
    predictions <- rbind(predictions, predictions_one_tumor_type)
    # }
    
  }
  # if(pan == "pan"){
  #   return(predictions_pan_cancer)
  # } else {
  return(predictions)
  # }
  
}

create_ROC_comparison <- function(tumor_type, model, feature_selection){
  predictions <- create_table_for_model("TP53_many_tumor_types", model, feature_selection, "NA")
  predictions <- predictions[predictions$Set == "test" & predictions$tumor_type == tumor_type,]
  pred <- prediction(predictions$pred_score, predictions$label)
  perf <- performance(pred, "tpr", "fpr")
  return(perf)
}

# create and return organized table of AUC values for every model/tumor type, ordered by number of labels.
Create_test_AUC_table_for_models_and_tumor_types = function(){
  # create tables for models 
  GCN_predictions <- create_table_for_model("TP53_many_tumor_types", "GCN", "NA", "NA")
  ELR_mutsigdb_predictions <- create_table_for_model("TP53_many_tumor_types", "ELR", "mutsigdb", "NA")
  RF_mutsigdb_predictions <- create_table_for_model("TP53_many_tumor_types", "RF", "mutsigdb", "NA")
  ELR_deseq_predictions <- create_table_for_model("TP53_many_tumor_types", "ELR", "deseq", "NA")
  RF_deseq_predictions <- create_table_for_model("TP53_many_tumor_types", "RF", "deseq", "NA")
  
  # create table of test set results for each tumor/model
  models <- c("GCN", "ELR_mutsigdb", "ELR_deseq", "RF_mutsigdb", "RF_deseq")
  tumor_types <- unique(GCN_predictions$tumor_type)
  
  test_AUC_models_tumor_types <- expand.grid(models,tumor_types)
  colnames(test_AUC_models_tumor_types) <- c("model", "tumor_type")
  test_AUC_models_tumor_types$test_AUC <- NA
  test_AUC_models_tumor_types$test_ACC <- NA
  test_AUC_models_tumor_types$test_prec <- NA
  test_AUC_models_tumor_types$test_recall <- NA
  test_AUC_models_tumor_types$val_AUC <- NA
  test_AUC_models_tumor_types$val_ACC <- NA
  test_AUC_models_tumor_types$val_prec <- NA
  test_AUC_models_tumor_types$val_recall <- NA
  test_AUC_models_tumor_types$overall_AUC <- NA
  test_AUC_models_tumor_types$overall_ACC <- NA
  test_AUC_models_tumor_types$n_samples <- NA
  test_AUC_models_tumor_types$percent_mutated <- NA
  test_AUC_models_tumor_types$log_n_samples <- NA
  
  for(model_full in models){
    
    feature_selection <- tolower(strsplit(model_full, "_")[[1]][2])
    model <- strsplit(model_full, "_")[[1]][1]
    
    for(tumor_type in tumor_types){
      
      temp_df <- get(paste0(model_full, "_predictions"))
      temp_df <- temp_df[temp_df$tumor_type == tumor_type,]
      model_tumor_type_bool_term <- test_AUC_models_tumor_types$tumor_type == tumor_type & 
        test_AUC_models_tumor_types$model == model_full
      
      test_ROC <- roc(predictor=temp_df$pred_score[temp_df$Set == "test"],
                      response=temp_df$label[temp_df$Set == "test"],
                      levels=rev(levels(as.factor(temp_df$label))))
      
      test_ACC <- sum(temp_df$label[temp_df$Set == "test"] == temp_df$prediction[temp_df$Set == "test"]) / 
        sum(temp_df$Set == "test")
      
      TP <- sum(temp_df$Set == "test" & temp_df$label == 1 & temp_df$prediction == 1)
      FP <- sum(temp_df$Set == "test" & temp_df$label == 0 & temp_df$prediction == 1)
      TN <- sum(temp_df$Set == "test" & temp_df$label == 0 & temp_df$prediction == 0)
      FN <- sum(temp_df$Set == "test" & temp_df$label == 1 & temp_df$prediction == 0)
      test_prec <- TP/(TP+FP)
      test_recall <- TP/(TP+FN)
      test_F1 <- 2*((test_prec*test_recall)/(test_prec+test_recall))
      # Matthews Correlation Coefficient
      test_mcc <- (TP*TN - FP*TN)/
        ( sqrt(
          (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)
        ) )
      test_F1[is.na(test_F1)] <- 0
      test_mcc[is.na(test_mcc)] <- -1
      
      
      val_ROC <-  roc(predictor=temp_df$pred_score[temp_df$Set == "val"],
                      response=temp_df$label[temp_df$Set == "val"],
                      levels=rev(levels(as.factor(temp_df$label))))
      
      val_ACC <- sum(temp_df$label[temp_df$Set == "val"] == temp_df$prediction[temp_df$Set == "val"]) / 
        sum(temp_df$Set == "val")
      
      TP <- sum(temp_df$Set == "val" & temp_df$label == 1 & temp_df$prediction == 1)
      FP <- sum(temp_df$Set == "val" & temp_df$label == 0 & temp_df$prediction == 1)
      TN <- sum(temp_df$Set == "val" & temp_df$label == 0 & temp_df$prediction == 0)
      FN <- sum(temp_df$Set == "val" & temp_df$label == 1 & temp_df$prediction == 0)
      val_prec <- TP/(TP+FP)
      val_recall <- TP/(TP+FN)
      val_F1 <- 2*((val_prec*val_recall)/(val_prec+val_recall))
      val_mcc <- (TP*TN - FP*TN)/
        ( sqrt(
          (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)
        ) )
      val_F1[is.na(val_F1)] <- 0
      val_mcc[is.na(val_mcc)] <- -1
      overall_ROC <-  roc(predictor=temp_df$pred_score,
                          response=temp_df$label,
                          levels=rev(levels(as.factor(temp_df$label))))
      
      overall_ACC <- sum(temp_df$label == temp_df$prediction) / 
        nrow(temp_df)
      
      
      test_AUC_models_tumor_types[model_tumor_type_bool_term,"test_AUC"] <- test_ROC$auc
      test_AUC_models_tumor_types[model_tumor_type_bool_term,"test_ACC"] <- test_ACC
      test_AUC_models_tumor_types[model_tumor_type_bool_term,"val_AUC"] <- val_ROC$auc
      test_AUC_models_tumor_types[model_tumor_type_bool_term,"val_ACC"] <- val_ACC
      test_AUC_models_tumor_types[model_tumor_type_bool_term,"overall_AUC"] <- overall_ROC$auc
      test_AUC_models_tumor_types[model_tumor_type_bool_term,"overall_ACC"] <- overall_ACC
      test_AUC_models_tumor_types[model_tumor_type_bool_term,"n_samples"] <- nrow(temp_df)
      test_AUC_models_tumor_types[model_tumor_type_bool_term,"log_n_samples"] <- log(nrow(temp_df))
      test_AUC_models_tumor_types[model_tumor_type_bool_term,"percent_mutated"] <- sum(temp_df$label == 1) / 
        nrow(temp_df)
      test_AUC_models_tumor_types[model_tumor_type_bool_term,"test_prec"] <- test_prec
      test_AUC_models_tumor_types[model_tumor_type_bool_term,"test_recall"] <- test_recall
      test_AUC_models_tumor_types[model_tumor_type_bool_term,"test_F1"] <- test_F1
      test_AUC_models_tumor_types[model_tumor_type_bool_term,"test_mcc"] <- test_mcc
      test_AUC_models_tumor_types[model_tumor_type_bool_term,"val_prec"] <- val_prec
      test_AUC_models_tumor_types[model_tumor_type_bool_term,"val_recall"] <- val_recall
      test_AUC_models_tumor_types[model_tumor_type_bool_term,"val_F1"] <- val_F1
      test_AUC_models_tumor_types[model_tumor_type_bool_term,"val_mcc"] <- val_mcc
    }
  }
  
  # exclude tumor types by downsampling criteria
  exclude_tumor_types = c("UCS", "KICH", "ACC","READ", "LAML", "ESCA")
  test_AUC_models_tumor_types <- test_AUC_models_tumor_types[!test_AUC_models_tumor_types$tumor_type %in% 
                                                               exclude_tumor_types,]
  
  # order levels of tumor_types by n_sample
  test_AUC_models_tumor_types$tumor_type <- as.character(test_AUC_models_tumor_types$tumor_type)
  test_AUC_models_tumor_types$tumor_type <- factor(test_AUC_models_tumor_types$tumor_type, 
                                                   levels =  unique(
                                                     test_AUC_models_tumor_types$tumor_type
                                                     [order(test_AUC_models_tumor_types$n_samples)]
                                                   ))
  test_AUC_models_tumor_types$model <- factor(test_AUC_models_tumor_types$model, levels = models)
  test_AUC_models_tumor_types <- test_AUC_models_tumor_types[order(rev(test_AUC_models_tumor_types$model)),]
  return(test_AUC_models_tumor_types)
}



### Data ###

# see colors. change name to: Set1/Set2/Set3
# look here for more: http://www.sthda.com/english/wiki/colors-in-r#using-rcolorbrewer-palettes
display.brewer.pal(n = 5, name = 'PuBu')
brw_dark2 <- brewer.pal(n = 3, name = "Dark2")

### 1. get AUC values from best_models.
best_models <- read.csv("//home//gil//Desktop//GCN//GCN_analysis//3.Models_comparison_grid//best_models_GCN_ELR_RF_multiple_cancer_types_1.2.csv", sep = ",", row.names = 1)
best_models <- best_models[!(best_models$model == "GCN"),]
# Integrate new GCN results
best_models_GCN <- read.csv("//home//gil//Desktop//GCN//GCN_analysis//3.Models_comparison_grid//best_models_TP53_GCN_multiple_cancer_types_3.0.csv", sep = ",")
best_models = rbind(best_models, best_models_GCN)
rm(best_models_GCN)

# add gene selection (deseq, mutsigdb) into models representation
best_models$model[best_models$model == "ELR" & best_models$gene_selection == "deseq"] <-"ELR_deseq" 
best_models$model[best_models$model == "ELR" & best_models$gene_selection == "mutsigdb"] <-"ELR_mutsigdb"
best_models$model[best_models$model == "LR" & best_models$gene_selection == "deseq"] <-"LR_deseq"
best_models$model[best_models$model == "LR" & best_models$gene_selection == "mutsigdb"] <-"LR_mutsigdb"
best_models$model[best_models$model == "RF" & best_models$gene_selection == "deseq"] <-"RF_deseq"
best_models$model[best_models$model == "RF" & best_models$gene_selection == "mutsigdb"] <-"RF_mutsigdb"
best_models$model <- factor(
  best_models$model, levels = c("GCN", "MLP", "ELR_mutsigdb", "ELR_deseq",
                                "LR_mutsigdb", "LR_deseq", "RF_mutsigdb", "RF_deseq"))

### ROC curves for Figure 2 ###

PATH_plots <- "//home//gil//Desktop//GCN//GCN_analysis//Figures_for_paper//Figure_2_models_comparison_metrices///"

for(tumor_type in c("BRCA", "LGG", "LUAD", "pan_cancer")){
  assign(paste0("perf_", tumor_type, "_GCN"), create_ROC_comparison(tumor_type, "GCN", "NA"))
  assign(paste0("perf_", tumor_type, "_ELR_mutsigdb"), create_ROC_comparison(tumor_type, "ELR", "mutsigdb"))
  assign(paste0("perf_", tumor_type, "_ELR_deseq"), create_ROC_comparison(tumor_type, "ELR", "deseq"))
  assign(paste0("perf_", tumor_type, "_RF_mutsigdb"), create_ROC_comparison(tumor_type, "RF", "mutsigdb"))
  assign(paste0("perf_", tumor_type, "_RF_deseq"), create_ROC_comparison(tumor_type, "RF", "deseq"))
  
  tiff(paste0(PATH_plots, tumor_type, "_ROC_3_Models_qualitative_3.0.tiff"),
       units="in", width=8, height=8, res=300)
  par(cex.axis=2, cex.lab = 2, mar=c(5,6,4,1))
  print(plot(get(paste0("perf_", tumor_type, "_GCN")), 
             col = brw_dark2[1], lwd = 7)) 
  # print(plot(get(paste0("perf_", tumor_type, "_ELR_deseq")), add = TRUE,
  #      col = brw_dark2[1], lwd = 7))
  # print(plot(get(paste0("perf_", tumor_type, "_RF_mutsigdb")), add = TRUE,
  #      col = brw_dark2[4], lwd = 7))
  print(plot(get(paste0("perf_", tumor_type, "_ELR_mutsigdb")), add = TRUE,
             col = brw_dark2[2], lwd = 7))
  print(plot(get(paste0("perf_", tumor_type, "_RF_deseq")), add = TRUE,
             col = brw_dark2[3], lwd = 7))
  GCN_legend <- paste0("GCN (AUC = ", round(100*best_models$test_auc[best_models$model == "GCN" &
                                                             best_models$cancer.type == tumor_type],2), "%)")
  ELR_msig_legend <- paste0("ELR MutSigDB (AUC = ", round(100*best_models$test_auc[best_models$model == "ELR_mutsigdb" &
                                                                           best_models$cancer.type == tumor_type],2), "%)")
  ELR_deseq_legend <- paste0("ELR Deseq (AUC = ", round(100*best_models$test_auc[best_models$model == "ELR_deseq" &
                                                                         best_models$cancer.type == tumor_type],2), "%)")
  RF_msig_legend <- paste0("RF MutSigDB (AUC = ", round(100*best_models$test_auc[best_models$model == "RF_mutsigdb" &
                                                                         best_models$cancer.type == tumor_type],2), "%)")
  RF_deseq_legend <- paste0("RF Deseq (AUC = ", round(100*best_models$test_auc[best_models$model == "RF_deseq" &
                                                                       best_models$cancer.type == tumor_type],2), "%)")
  
  legend("bottomright",
         legend=c(GCN_legend, ELR_msig_legend, RF_deseq_legend),
         col=brw_dark2, 
         lwd=7, cex =1.5, xpd = TRUE, horiz = F)
  
  
  
  dev.off()
  
}



### BARPLOT OF DIFFERENT METRICES for Figure 2 ###
# create organized table for AUC in all tumor-type/model pairs.
test_AUC_models_tumor_types = Create_test_AUC_table_for_models_and_tumor_types()


best_models_bar_plot <- test_AUC_models_tumor_types[,c("model", "tumor_type", "test_AUC", 
                                                       "test_ACC", "test_F1", "n_samples", "percent_mutated")]
best_models_bar_plot <- pivot_longer(best_models_bar_plot, cols=c("test_AUC", "test_ACC", "test_F1"), names_to = "metric", values_to = "metric_value")
best_models_bar_plot <- best_models_bar_plot[best_models_bar_plot$tumor_type %in% c("pan_cancer", "BRCA", "LGG", "LUAD"),]
best_models_bar_plot <- best_models_bar_plot[best_models_bar_plot$model %in% c("GCN", "ELR_mutsigdb", "RF_deseq"),]
best_models_bar_plot$metric_value <- 100*round(best_models_bar_plot$metric_value, 4)
best_models_bar_plot$percent_mutated <- 100*round(best_models_bar_plot$percent_mutated, 4)


best_models_bar_plot$tumor_type <- factor(best_models_bar_plot$tumor_type, levels=c( "pan_cancer", "BRCA", "LGG", "LUAD"))
levels(best_models_bar_plot$tumor_type) <- c("Pan Cancer", "BRCA", "LGG", "LUAD")

tiff(paste0("//home//gil//Desktop//GCN//GCN_analysis//Figures_for_paper//Figure_2_models_comparison_metrices///Barplot_by_tumor_type_AUC_Acc_F1_3.0.tiff"), units="in", width=9, height=8, res=300)
ggplot(best_models_bar_plot, aes(x=metric, y=metric_value, fill=model)) +
  geom_bar(stat="identity",  width=0.9, position=position_dodge(0.9)) + 
  scale_fill_manual(labels = c("GCN", "ELR MutSigDB", "RF Deseq"), values = brw_dark2) +theme_bw()+
   theme(legend.title=element_text(size=18))+ theme(legend.text=element_text(size=18))+
  geom_text(aes(label=metric_value), hjust = 1.2, angle = 90, position = position_dodge(0.9), color="white", size=6) + 
  labs(x = "Metric", y = "Value") + facet_wrap(~tumor_type, labeller = labeller("pan_cancer" = "Pan Cancer", "LUAD" = "LUAD",
                                                                                "LGG" = "LGG", "BRCA" = "BRCA"))+
  scale_x_discrete(labels = c("Test ACC", "Test AUC", "Test F1")) +font("y.text", size = 18) +font("x.text", size = 15) +
  font("xy.title", size = 19) + guides(fill=guide_legend(title="Model")) + theme(strip.text.x = element_text(size=20))
  
dev.off()



### All tumor types models compare by F1 / ACC/ AUC for Figure 3 ###
set = "val"
# try only 3 models. possible with all 5.
# select the variable: score
# score = "F1"
score = "AUC"

test_AUC_models_tumor_types <- test_AUC_models_tumor_types[test_AUC_models_tumor_types$model %in% c("GCN", "ELR_mutsigdb", "RF_deseq"),]
tiff(paste0("//home//gil//Desktop//GCN//GCN_analysis//Figures_for_paper//Figure_3_Downsampling_and_all_tumor_types///Models_comparison_by_", 
            set, "_", score, "_as.numeric_factor_tumor_types_3.0.tiff"), units="in", width=24, height=8, res=300)
p <- ggplot(data = test_AUC_models_tumor_types,
            aes(x = as.numeric(tumor_type), y = get(paste0(set, "_", score)),  colour = model, shape = model)) +
  geom_point(size = 10) + labs(x = "Tumor Type", 
                               y = paste0("Test ", score)) + theme_bw()+
  stat_smooth(method = "loess", se = F, size = 5) + scale_color_manual(values = brw_dark2)+
  theme(plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  font("x.text", size = 30) +font("y.text", size = 40) +
  font("legend.title", size = 40) + font("legend.text", size = 40) + 
  font("xy.title", size = 40) #+ labs(color=l, shape=l)#+ guides(col=guide_legend(title="Models")) #+
# scale_x_discrete(labels=unique(test_AUC_models_tumor_types$tumor_type))
print(p)
dev.off()


### Log Sample Size Barplot for Figure 3 ###

samples_barplot <- test_AUC_models_tumor_types[,c("tumor_type", "n_samples")]
samples_barplot <- samples_barplot[!duplicated(samples_barplot),]
samples_barplot <- samples_barplot[order(samples_barplot$n_samples),]
levels(samples_barplot$tumor_type)[24] <- "Pan Cancer"
samples_barplot$tumor_type[24] <- "Pan Cancer"

# change green to blue?
tiff(paste0("//home//gil//Desktop//GCN//GCN_analysis//Figures_for_paper//Figure_3_Downsampling_and_all_tumor_types///Barplot_num_samples_tumor_types.tiff"), units="in", width=16, height=8, res=300)
ggplot(samples_barplot, aes(x=tumor_type, y=log(n_samples), fill=log(n_samples))) + 
  geom_bar(stat="identity",  width=0.95, position=position_dodge(0.9)) + 
  scale_fill_distiller(type = "div", palette = "",direction = 1) + theme_void() + font("x.text", , size = 35, angle = 90) + 
  labs(x = "Tumor Type", y = "Sample Size") + font("x.title", size = 15) + font("y.title", size = 35, angle = 90)

dev.off()

