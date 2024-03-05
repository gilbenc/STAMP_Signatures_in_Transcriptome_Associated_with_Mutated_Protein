### 
### Survival analysis for TCGA expression project
###
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(survminer)
library(survival)
library("RColorBrewer")
library(ggfortify)

PATH_datasets = "//home//gil//Desktop//GCN//datasets//Data_downloaded_directly_from_source//"
PATH_pred_tables = "//home//gil//Desktop//GCN//GCN_analysis//1.Prediction_tables//"
PATH_survival_figures = "//home//gil//Desktop//GCN//GCN_analysis//5.Linear_Score_and_Survival//"



### Functions ###

# create a function to calculate the logit of a value
logit <- function(x) {
        logit_value <- log(x/(1-x))
        
        return(logit_value)
}

# Simple vector normalization (0 to 1) based on min and max values.
normalize <- function(x) { return( (x- min(x)) / (max(x) - min(x)) ) }

# get vector of pred_scores and return linear score.
create_linear_score <- function(x, model) {
  if(model != "GCN"){
    x <- normalize(x)
  }
  x <- sapply(x, logit)
  # Provide limits to logit when equals infinity
  # Keep as larger than largest value for comparison
  x[x == Inf] <- max(x[x != Inf]) + 2
  x[x == -Inf] <- min(x[x != -Inf]) - 2
  
  return(x)
}

# input: pred_table with columns: "Sample ID", "label", "prediction", "linear_score"
# output: combined df with months and status for survival
df_for_survival <- function(pred_table) {
  
  ### load TCGA Clinical data ###
  clinical.path = file.path(paste0(PATH_datasets, "Cbioportal_Pan_Cacner_WmutProfile_clinical_data.tsv"), fsep = .Platform$file.sep)
  clinical_data<- read.csv(clinical.path, sep="\t", fileEncoding="UTF-8-BOM", header=T, stringsAsFactors = F, check.names = F)
  if(tumor_type == "pan_cancer"){
    survival_data <- clinical_data[,c('Sample ID', 'Overall Survival (Months)', 'Overall Survival Status')]
    # survival_data <- clinical_data[,c('Sample ID', 'Overall Survival (Months)', 'Overall Survival Status', "Diagnosis Age", "TCGA PanCanAtlas Cancer Type Acronym","Mutation Count")]
  } else {
    survival_data <- clinical_data[clinical_data$`TCGA PanCanAtlas Cancer Type Acronym` == tumor_type,c('Sample ID', 'Overall Survival (Months)', 'Overall Survival Status')]
    # survival_data <- clinical_data[clinical_data$`TCGA PanCanAtlas Cancer Type Acronym` == tumor_type,c('Sample ID', 'Overall Survival (Months)', 'Overall Survival Status', "Diagnosis Age", "TCGA PanCanAtlas Cancer Type Acronym","Mutation Count")]
  }
  # merge with models
  pred_table_survival <- merge(pred_table, survival_data, 'Sample ID')
  
  # survival analysis
  months <- pred_table_survival$`Overall Survival (Months)`
  status <- pred_table_survival$`Overall Survival Status`
  
  label <- pred_table_survival$label
  binary_pred <- pred_table_survival$prediction
  linear_pred <- pred_table_survival$linear_score
  
  # tumor type column
  # tumor_type_col <- pred_table_survival$`TCGA PanCanAtlas Cancer Type Acronym`
  # age <- pred_table_survival$`Diagnosis Age`
  # Mutation_count <- pred_table_survival$`Mutation Count`
  
  # 
  # MSI_sensor <-pred_table_survival$`MSIsensor Score`
  # MSI_mantis <- pred_table_survival$`MSI MANTIS Score`
  # categorical_pred <- # to be modified
  
  status <- as.integer(status == "1:DECEASED")
  dataframe <- data.frame(months, status, binary_pred, linear_pred, label)
  # dataframe <- data.frame(months, status, binary_pred, linear_pred, label, tumor_type_col, age, Mutation_count)
  # quantiles_label <- as.numeric(dataframe$linear_pred > quantile(dataframe$linear_pred, 0.25))
  # quantiles_label[dataframe$linear_pred > quantile(dataframe$linear_pred, 0.5)] <- 2
  # quantiles_label[dataframe$linear_pred > quantile(dataframe$linear_pred, 0.75)] <- 3
  # dataframe$quantiles_label <- quantiles_label
  
  return(dataframe)
}

# input: dataframe with columns: "Sample ID", "label", "prediction", 
# "linear_score", "months" and "status".
# output: the "row" for the p-values dataframe
create_survival_pvalues_row <- function(dataframe){
  
  # label
  fit_coxph <- coxph(Surv(months, status) ~ label, data = dataframe)
  sum_fit_coxph <- summary(fit_coxph)
  label_logtest_p <- sum_fit_coxph$logtest["pvalue"][[1]]
  # pred
  fit_coxph <- coxph(Surv(months, status) ~ binary_pred, data = dataframe)
  sum_fit_coxph <- summary(fit_coxph)
  binary_logtest_p <- sum_fit_coxph$logtest["pvalue"][[1]]
  # linear
  fit_coxph <- coxph(Surv(months, status) ~ linear_pred, data = dataframe)
  sum_fit_coxph <- summary(fit_coxph)
  linear_logtest_p <- sum_fit_coxph$logtest["pvalue"][[1]]
  
  # pred quantiles
  up_q <- quantile(dataframe$linear_pred, 0.75)
  down_q <- quantile(dataframe$linear_pred, 0.25)
  dataframe_quantile <- dataframe[(dataframe$linear_pred > up_q) | (dataframe$linear_pred < down_q),]
  if(nrow(dataframe_quantile) == 0) { 
    quantile_logtest_p <- 0
  } else {
    fit_coxph <- coxph(Surv(months, status) ~ binary_pred, data = dataframe_quantile)
    sum_fit_coxph <- summary(fit_coxph)
    quantile_logtest_p <- sum_fit_coxph$logtest["pvalue"][[1]] 
  }
  
  is_pred_better <- (linear_logtest_p < label_logtest_p)
  
  survival_pvalues_row <- c(label_logtest_p, binary_logtest_p,
                            linear_logtest_p, quantile_logtest_p, is_pred_better)
  return(survival_pvalues_row)
  
}

# input: dataframe with linear score and survival values.
# output: optimal cutoff by survival distinction p-val.
optimize_cutoff_by_coxph <- function(dataframe){
  best_cutoff <- 0
  best_p <- 1
  for(i in 1:nrow(dataframe)){
    cutoff <- dataframe$linear_pred[i]
    dataframe$binary_pred <- as.numeric(dataframe$linear_pred > cutoff)
    fit_coxph <- coxph(Surv(months, status) ~ binary_pred, data = dataframe)
    sum_fit_coxph <- summary(fit_coxph)
    logtest_p <- sum_fit_coxph$logtest["pvalue"][[1]]  
    if(logtest_p < best_p){
      best_cutoff <- cutoff
      best_p <- logtest_p
    }
  }
  return(best_cutoff)
}

# input: dataframe with linear score a
# output: optimal cutoff by accuracy (compared to label). Without!! using test.
optimize_cutoff_by_label <- function(pred_table){
  pred_table <- pred_table[pred_table$Set != "test",]  
  best_cutoff <- 0
  best_acc <- 0
  for(i in 1:nrow(pred_table)){
    cutoff <- pred_table$linear_score[i]
    pred_table$prediction <- as.numeric(pred_table$linear_score > cutoff)
    acc <- sum(pred_table$prediction == pred_table$label)
    
    if(acc > best_acc){
      best_cutoff <- cutoff
      best_acc <- acc
    }
  }
  return(best_cutoff)
}

# dataframe = df required for survival
# label_type = label / pred / quantiles. (NOT LINEAR, because it can't show binary curves). 
# for quantiles, show 4 groups by .25, .5, .75 separation
create_survival_curve <- function(dataframe, data_type, gene, tumor_type, model, 
                                  feature_selection, quantiles, label_type, cutoff="pred") {
  if(label_type == "quantiles"){
    quantiles_label <- as.numeric(dataframe$linear_pred > quantile(dataframe$linear_pred, 0.25))
    quantiles_label[dataframe$linear_pred > quantile(dataframe$linear_pred, 0.5)] <- 2
    quantiles_label[dataframe$linear_pred > quantile(dataframe$linear_pred, 0.75)] <- 3
    dataframe$quantiles_label <- quantiles_label
    
    
    fit <- survfit(Surv(months, status) ~ quantiles_label, data = dataframe)
    print(fit)
    coxph(Surv(months, status) ~ quantiles_label, data = dataframe)
    
    fit_coxph_0_1 <- coxph(Surv(months, status) ~ quantiles_label, data = dataframe[dataframe$quantiles_label == 0 | dataframe$quantiles_label == 1,])
    summary(fit_coxph_0_1)
    
    fit_coxph_1_2 <- coxph(Surv(months, status) ~ quantiles_label, data = dataframe[dataframe$quantiles_label == 1 | dataframe$quantiles_label == 2,])
    summary(fit_coxph_1_2)
    
    fit_coxph_2_3 <- coxph(Surv(months, status) ~ quantiles_label, data = dataframe[dataframe$quantiles_label == 2 | dataframe$quantiles_label == 3,])
    summary(fit_coxph_2_3)
    

    
  }
  if(label_type == "label"){
    fit <- survfit(Surv(months, status) ~ label, data = dataframe)
    fit_coxph <- coxph(Surv(months, status) ~ label, data = dataframe)
    sum_fit_coxph <- summary(fit_coxph)
    logtest_p <- sum_fit_coxph$logtest["pvalue"][[1]]  
  }
  if(label_type == "pred"){
    if(quantiles){
      up_q <- quantile(dataframe$linear_pred, 0.75)
      down_q <- quantile(dataframe$linear_pred, 0.25)
      dataframe <- dataframe[dataframe$linear_pred > up_q | dataframe$linear_pred < down_q,]
      if(nrow(dataframe) == 0) return()
    }
    fit <- survfit(Surv(months, status) ~ binary_pred, data = dataframe)
    fit_coxph <- coxph(Surv(months, status) ~ binary_pred, data = dataframe)
    sum_fit_coxph <- summary(fit_coxph)
    logtest_p <- sum_fit_coxph$logtest["pvalue"][[1]]  
  }
  if(label_type == "linear"){
    fit <- survfit(Surv(months, status) ~ linear_pred, data = dataframe)
    fit_coxph <- coxph(Surv(months, status) ~ linear_pred, data = dataframe)
    sum_fit_coxph <- summary(fit_coxph)
    logtest_p <- sum_fit_coxph$logtest["pvalue"][[1]]
    
  }
  
  
  
  if(length(fit$n) != 2 && label_type != "quantiles") return()
  
  if(label_type == "pred") {
    legend_labs = c("ND", "D")
    legend_title = "Prediction" 
    } else {
      legend_labs = c("NM", "M")
      legend_title = "Label"
    }
  if(label_type == "quantiles") {
    legend_labs = c("Quant 1", "Quant 2", "Quant 3", "Quant 4")
    legend_title = "Quantiles"
  }
  title = paste0("Survival curve ", data_type, " ", gene, " ", tumor_type, " ", model)
  
  if(quantiles){
    subtitle = paste0(feature_selection, " quantiles: ", quantiles, "\nlog test p: ", logtest_p)  
  } else {
    if(label_type == "quantiles"){
      subtitle = ""
      palette_ggsurv = c("#2c7bb6", "#abd9e9","#fdae61" , "#d7191c")
    } else subtitle = paste0(feature_selection, " ", label_type, "\nlog test p: ", logtest_p)  
  }
  if(label_type != "quantiles"){
    logtest_p <- format(logtest_p, digits = 5)  
    palette_ggsurv = c("#67a9cf", "#ef8a62")# "#377EB8", "#F781BF"),
  }
  
  # create dir for this tumor type
  # dir.create(paste0(PATH_survival_figures, data_type, "//", tumor_type), showWarnings = F) 
  PATH_save_survival_curve <- paste0(PATH_survival_figures, data_type, "//",model,"//", tumor_type)
  dir.create(file.path(PATH_survival_figures, paste0(data_type, "//",model,"//")), showWarnings = FALSE)
  dir.create(file.path(PATH_survival_figures, paste0(data_type, "//",model,"//",tumor_type,"//")), showWarnings = FALSE)
  
  # calculate x axis limit
  x_axis_limit = 100+max(fit$time) - max(fit$time)%%100
  
  
  
  
  ggsurv <- ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    # title = title,
    # subtitle = subtitle,
    size = 10,
    data = dataframe,             # data used to fit survival curves.
    # risk.table = TRUE,       # show risk table.
    # pval = TRUE,             # show p-value of log-rank test.
    conf.int = F,         # show confidence intervals for
    # point estimates of survival curves.
    palette = palette_ggsurv,# "#377EB8", "#F781BF"),
    xlim = c(0,x_axis_limit),         # present narrower X axis, but not affect
    # survival estimates.
    xlab = "Time in months",# customize X axis label.
    break.time.by = 50,     # break X axis in time intervals by 500.
    ggtheme = theme_light(), # customize plot and risk table with a theme.
    # risk.table.y.text.col = T,# colour risk table text annotations.
    # risk.table.height = 0.25, # the height of the risk table
    # risk.table.y.text = FALSE,# show bars instead of names in text annotations
    # in legend of risk t)able.
    # ncensor.plot = TRUE,      # plot the number of censored subjects at time t
    # ncensor.plot.height = 0.25,
    # conf.int.style = "step",  # customize style of confidence intervals
    
    # surv.median.line = "hv",  # add the median survival pointer.
    legend.labs = legend_labs
    ,    # change legend labels.
    legend = "right"
  ) 
  if(label_type != "quantiles") ggsurv$plot <- ggsurv$plot + ggplot2::annotate("text", 
                                                    x = x_axis_limit-50, y = 0.9, # x and y coordinates of the text
                                                    label = paste0("Log Test \np = ", logtest_p), size = 15)
  ggsurv <- ggpar(
    ggsurv,
    # font.title    = c(16, "bold", "darkblue"),
    # font.subtitle = c(15, "bold.italic", "purple"),
    # font.caption  = c(30, "#E41A1C", "#4DAF4A", "#377EB8"),
    font.x        = c(40, "plain", "black"),
    font.y        = c(40, "plain", "black"),
    font.xtickslab = c(40, "plain", "black"),
    font.ytickslab = c(40, "plain", "black"),
    font.legend =  c(40, "plain", "black"),
    
    legend.title = legend_title
  )
  
  if(quantiles) {
    tiff(paste0(PATH_save_survival_curve, "//Survival_", label_type, "_", gene, "_",  tumor_type, "_",
                model, feature_selection, "_quantiles_logtest_pval_", logtest_p, ".tiff"), units="in", width=15, height=8, res=300)
    print(ggsurv)
    dev.off()
  } else { 
    if(label_type == "quantiles") {
      tiff(paste0(PATH_save_survival_curve, "//Survival_", label_type, "_", gene, "_", tumor_type, "_",
                  model, feature_selection, "4_quantiles.tiff"), units="in", width=15, height=8, res=300)
      print(ggsurv)
      dev.off()
    } else {
        tiff(paste0(PATH_save_survival_curve, "//Survival_", label_type, "_", gene, "_", tumor_type, "_",
                    model, feature_selection, "_", cutoff,"_logtest_pval_", logtest_p, ".tiff"), units="in", width=15, height=8, res=300)
        print(ggsurv)
        dev.off()
    }
  }
}

# get relevant values and create a row and curves for survival analysis
create_survival_row_and_curve <- function(model, data_type, tumor_type, gene, create_label_curve) {
  # create path and reduce files to the one I want
  PATH <- paste0(PATH_pred_tables, data_type, "//")
  files <- list.files(paste0(PATH, model, "//"))
  files <- files[grepl(files, pattern = gene)]
  files <- files[grepl(files, pattern = tumor_type)]
  if(model == "ELR") {
    feature_selection <- "mutsigdb"
    if(data_type != "Pathways") files <- files[grepl(files, pattern = "mutsigdb")] 
  }
  if(model == "RF") {
    if(data_type == "Pathways") return(0)
    feature_selection <- "deseq"
    files <- files[grepl(files, pattern = "deseq")]
  }
  if(model == "GCN") feature_selection <- ""
  if(length(files) == 0) return(0)
  if(length(files) == 2) {
      if(gene == "NOTCH") {
        files <- files[!grepl(files, pattern = "NOTCH1")] 
      } else {
        files <- files[!grepl(files, pattern = "pathway")]
      }
    }
  
  # read and prepare the file
  pred_table <- read.csv(paste0(PATH, model, "//", files), row.names = "Sample.ID")
  pred_table <- pred_table[,-1]
  pred_table$"Sample ID" <- row.names(pred_table)
  
  pred_table$linear_score <- create_linear_score(x = pred_table$pred_score, model = model)
  
  # dataframe for survival analysis
  dataframe <- df_for_survival(pred_table)    
  

  
  # get survival p-values for p-values table
  survival_pvalues_row <- create_survival_pvalues_row(dataframe)
  survival_pvalues_row <- c(data_type, tumor_type, gene, model, feature_selection, survival_pvalues_row)
  

  
  # generate a survival curve
  if(create_label_curve) create_survival_curve(dataframe, data_type, gene, tumor_type, model, feature_selection, quantiles = F, label_type = "label")
  
  create_survival_curve(dataframe, data_type, gene, tumor_type, model, feature_selection, quantiles = F, label_type = "pred", cutoff="pred")
  
  # generate plot with coxph optimal cutoff
  # cutoff <- optimize_cutoff_by_coxph(dataframe)
  # dataframe$binary_pred <- as.numeric(dataframe$linear_pred > cutoff)
  # create_survival_curve(dataframe, data_type, gene, tumor_type, model, feature_selection, quantiles = F, label_type = "pred", cutoff="coxph")
  
  # create p values for cox 
  # fit_coxph <- coxph(Surv(months, status) ~ binary_pred, data = dataframe)
  # sum_fit_coxph <- summary(fit_coxph)
  # coxph_logtest_p <- sum_fit_coxph$logtest["pvalue"][[1]]
  # 
  
  # # generate plot with acc optimal accuraccy
  # cutoff <- optimize_cutoff_by_label(dataframe)
  # dataframe$binary_pred <- as.numeric(dataframe$linear_pred > cutoff)
  # create_survival_curve(dataframe, data_type, gene, tumor_type, model, feature_selection, quantiles = F, label_type = "pred", cutoff="label")
  # 
  # # create p values for acc
  # fit_coxph <- coxph(Surv(months, status) ~ binary_pred, data = dataframe)
  # sum_fit_coxph <- summary(fit_coxph)
  # acc_logtest_p <- sum_fit_coxph$logtest["pvalue"][[1]]
  # 
  # survival_pvalues_row = c(survival_pvalues_row, coxph_logtest_p)
  

  
  # create_survival_curve(dataframe, data_type, gene, tumor_type, model, feature_selection, quantiles = T, label_type = "pred")
  
  return(survival_pvalues_row)
  
}

# get list of relevant genes from best_models tables in pathway analysis
get_genes_list_from_pathway_analysis <- function(){
  pathway_path <- "//home//gil//Desktop//GCN//GCN_analysis//6.Cell_pathway_labeling_grid//"
  pathway_best_models <- read.csv(paste0(pathway_path, "best_models_ELR_GCN_Cell_paper_labels_PATHWAY.csv"))
  gene_best_models <- read.csv(paste0(pathway_path, "best_models_ELR_GCN_Cell_paper_labels_GENE.csv"))
  
  
  # adjust corrected GCN results
  pathway_best_models = pathway_best_models[!(pathway_best_models$model == "GCN"),]
  gene_best_models = gene_best_models[!(gene_best_models$model == "GCN"),]
  
  pathway_best_models_GCN <- read.csv(paste0(pathway_path, "best_models_GCN_Cell_paper_labels_PATHWAY_FIXED_ADJ.csv"))
  gene_best_models_GCN <- read.csv(paste0(pathway_path, "best_models_GCN_Cell_paper_labels_GENE_FIXED_ADJ.csv"))
  pathway_best_models = rbind(pathway_best_models, pathway_best_models_GCN)
  gene_best_models = rbind(gene_best_models, gene_best_models_GCN)
  rm(gene_best_models_GCN, pathway_best_models_GCN)
  
  pathway_best_models$gene[pathway_best_models$gene %in% c("MYC", "TP53")] <- paste0(pathway_best_models$gene[pathway_best_models$gene %in% c("MYC", "TP53")], "_pathway")
  genes <- unique(c("TP53", pathway_best_models$gene, gene_best_models$gene))
  return(genes)
}


##############
###  MAIN  ###
##############
clinical.path = file.path(paste0(PATH_datasets, "Cbioportal_Pan_Cacner_WmutProfile_clinical_data.tsv"), fsep = .Platform$file.sep)
clinical_data<- read.csv(clinical.path, sep="\t", fileEncoding="UTF-8-BOM", header=T, stringsAsFactors = F, check.names = F)
tumor_types <- unique(clinical_data$`TCGA PanCanAtlas Cancer Type Acronym`)
tumor_types <- c(tumor_types, "pan_cancer")

# need_survival_pvals = F
keep_original_pred = T

### Code for creating curves to 6 tumor types ###
### Adjust cutoff, insert p-value into figure ###
for(tumor_type in c("HNSC", "LGG", "LUAD", "UCEC", "PAAD", "COAD")){
  PATH <- paste0(PATH_pred_tables, data_type, "//")
  files <- list.files(paste0(PATH, model, "//"))
  files <- files[grepl(files, pattern = gene)]
  files <- files[grepl(files, pattern = tumor_type)]
  if(model == "ELR") {
    feature_selection <- "mutsigdb"
    if(data_type != "Pathways") files <- files[grepl(files, pattern = "mutsigdb")] 
  }
  if(model == "RF") {
    if(data_type == "Pathways") return(0)
    feature_selection <- "deseq"
    files <- files[grepl(files, pattern = "deseq")]
  }
  if(model == "GCN") feature_selection <- ""
  if(length(files) == 0) return(0)
  if(length(files) == 2) files <- files[!grepl(files, pattern = "pathway")]
  
  # read and prepare the file
  pred_table <- read.csv(paste0(PATH, model, "//", files), row.names = "Sample.ID")
  pred_table <- pred_table[,-1]
  pred_table$"Sample ID" <- row.names(pred_table)
  
  pred_table$linear_score <- create_linear_score(x = pred_table$pred_score, model = model)
  
  ### Create cutoff by label 
  # pred_table$prediction_preliminary <- pred_table$prediction
  # cutoff <- optimize_cutoff_by_label(pred_table)
  # print(paste0(tumor_type, " cutoff = ", cutoff))
  if(tumor_type == "HNSC") pred_table$prediction = as.numeric(pred_table$linear_score > 0)
  if(!keep_original_pred)  pred_table$prediction <- as.numeric(pred_table$linear_score > cutoff)
  
  pred_table$binary_pred = pred_table$prediction
  # input: pred table for specific model
  # output: column for new prediction column based on cutoff optimization
  # optimize using TRAIN & VAL ONLY.
  
  # dataframe for survival analysis
  dataframe <- df_for_survival(pred_table)    
  
  # get survival p-values for p-values table
  survival_pvalues_row <- create_survival_pvalues_row(dataframe)
  survival_pvalues_row <- c(data_type, tumor_type, gene, model, feature_selection, survival_pvalues_row)
  
  
  
  # generate a survival curve
  # if(create_label_curve) create_survival_curve(dataframe, data_type, gene, tumor_type, model, feature_selection, quantiles = F, label_type = "label")
  # dataframe <- dataframe[order(dataframe$linear_pred, decreasing = F),]
  # dataframe$binary_pred[1:round(nrow(dataframe)/2)] <- 0
  # dataframe$binary_pred[((round(nrow(dataframe))/2)+1):nrow(dataframe)] <- 1
  
  # cutoff <- optimize_cutoff_by_coxph(dataframe)
  # dataframe$binary_pred <- as.numeric(dataframe$linear_pred > cutoff)
  # print(paste0(tumor_type, " cutoff = ", cutoff, " n(D) = ", sum(dataframe$binary_pred == 1), " n(ND) = ",
  #              sum(dataframe$binary_pred == 0)))
  # create_survival_curve(dataframe, data_type, gene, tumor_type, model, feature_selection, quantiles = F, label_type = "label")
  create_survival_curve(dataframe, data_type, gene, tumor_type, model, feature_selection, quantiles = F, label_type = "pred")
  # create_survival_curve(dataframe, data_type, gene, tumor_type, model, feature_selection, quantiles = T, label_type = "pred")
  
}
