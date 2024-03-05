library(patchwork)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

display.brewer.pal(n = 5, name = 'PuBu')
brw_dark2 <- brewer.pal(n = 3, name = "Dark2")
brw_set2 <- brewer.pal(n = 5, name = "Set2")


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
  
  x[x == Inf] <- max(x[x != Inf]) + 2
  x[x == -Inf] <- min(x[x != -Inf]) - 2
  
  return(x)
}



WES_SBS_sigs = read.csv("//home//gil//Desktop//GCN//datasets//Data_downloaded_directly_from_source//TCGA_WES_SBS_sigs_realtive.csv", 
                        row.names = 1, stringsAsFactors = F)

WES_SBS_sigs <- as.data.frame(t(WES_SBS_sigs), stringsAsFactors = FALSE)

mut_sig_TCGA_IDs <- row.names(WES_SBS_sigs)
mut_sig_TCGA_IDs <- substr(mut_sig_TCGA_IDs, 1, 15)
row.names(WES_SBS_sigs) <- gsub('\\.', "-", mut_sig_TCGA_IDs)
WES_SBS_sigs$'Sample.ID' <- row.names(WES_SBS_sigs)

### Create boxplot and regression plots of interesting cancer type/gene/sbs
PATH <- "//home//gil//Desktop//GCN//GCN_analysis//"
PATH_figures <- paste0(PATH, "Figures_for_paper//Figure_4_Survival")
PATH_figures_regression = paste0(PATH_figures, "//regression_plots//")
PATH_pred_tables = "//home//gil//Desktop//GCN//GCN_analysis//1.Prediction_tables//"
model = "GCN"
SBS_tumor_type_gene_list = list(  c("SBS44", "COAD", "APC"), c("SBS20", "STAD", "NOTCH"))

## Regression plots
colors1 = c("#e6ab02","#7570b3","#66a61e")
colors2=c("#FFE979", "#c7c2fc", "#bbe889")
df = as.data.frame(c("model", "SBS", "tumor_type", "gene"))
df= t(df)
colnames(df) = df[1,]


spearman_rho = c()
spearman_pval = c()
wilc_label_pval = c()
wilc_pred_pval = c()
models = c()
SBS_s = c()
tumor_types = c()
genes = c()

#### ADJUST COLORS! (I only inserted BRCA-SBS13 colors)
i=1
j=1
for(model in c("ELR")) {#, "ELR", "RF")){

  for(SBS_tumor_type_gene in SBS_tumor_type_gene_list){
    
    ### Load model

    SBS = SBS_tumor_type_gene[1]
    tumor_type = SBS_tumor_type_gene[2]
    gene = SBS_tumor_type_gene[3]
    if(gene == "TP53") {data_type = "TP53_many_tumor_types"} else {data_type = "Pathways"}
    if(model == "ELR") {feature_selection = "mutsigdb"}
    if(model == "RF") {feature_selection = "deseq"}
  
    files <- list.files(paste0(PATH_pred_tables, data_type,"//", model, "//"))
    files_tumor_type <- files[grepl(files, pattern = tumor_type)]
    if(gene != "TP53") {files_tumor_type_gene <- files_tumor_type[grepl(files_tumor_type, pattern = gene)]} else {files_tumor_type_gene <- files_tumor_type[grepl(files_tumor_type, pattern = "TP53")]}
    if(model != "GCN" & data_type != "Pathways") {files_tumor_type_gene <- files_tumor_type_gene[grepl(files_tumor_type_gene, pattern = feature_selection)]}
    if(length(files_tumor_type_gene) == 0) next
    model_preds <- read.csv(paste0(PATH_pred_tables, data_type,"//", model, "//", files_tumor_type_gene))
    model_preds <- model_preds[,-1]
    
    model_preds <- merge(model_preds, WES_SBS_sigs[,c(SBS, "Sample.ID")], by = "Sample.ID")
    model_preds$linear_score <- create_linear_score(model_preds$pred_score, model)
    model_preds[,SBS] <- as.numeric(model_preds[,SBS])
    model_preds$tumor_type = tumor_type
    
    
    
    col_SBS1 = colors1[1]
    col_SBS2 = colors2[1]
    if((SBS == "SBS2" & tumor_type == "HNSC" & gene == "TP53") | (SBS == "SBS44" & tumor_type == "COAD" & gene == "APC")){
      col_SBS1 = colors1[2]
      col_SBS2 = colors2[2]
    }
    if((SBS == "SBS3" & tumor_type == "BRCA" & gene == "TP53") | (SBS == "SBS20" & tumor_type == "STAD" & gene == "NOTCH")){
      col_SBS1 = colors1[3]
      col_SBS2 = colors2[3]
    }
    
    
    
    ### Spearman for linear score
    test = cor.test(x=model_preds$linear_score,y = model_preds[,SBS], method="pearson")
    print(test)
    test = cor.test(x=model_preds$linear_score,y = model_preds[,SBS], method="spearman")
    print(test)
    spearman_rho = c(spearman_rho, test$estimate[[1]])
    spearman_pval = c(spearman_pval, test$p.value[[1]])
    
    ### wilcoxon for label and for pred
    wilc_test_label = wilcox.test(model_preds[,SBS]~model_preds$label)
    wilc_test_pred = wilcox.test(model_preds[,SBS]~model_preds$prediction)
    
    wilc_label_pval = c(wilc_label_pval, wilc_test_label$p.value)
    wilc_pred_pval = c(wilc_pred_pval, wilc_test_pred$p.value)
    models = c(models, model)
    SBS_s = c(SBS_s, SBS)
    tumor_types = c(tumor_types, tumor_type)
    genes = c(genes, gene)
    
    ## Regression plot
    
    const = 300
    if((tumor_type == "BRCA" & gene == "TP53" & SBS %in% c("SBS13", "SBS2")) |
       (SBS == "SBS44" & tumor_type == "COAD" & gene == "APC")) const = 500
    if(SBS == "SBS20" & tumor_type == "STAD" & gene == "NOTCH") const = 1000
    bins_const = 100
    if(tumor_type == "pan_cancer"){
      const = 200
      bins_const = 2000
      if(SBS %in% c("SBS2", "SBS13")) const = 500
      tumor_type = "Pan cancer"
    }
    if(tumor_type == "HNSC" & gene == "TP53" & SBS == "SBS2"){
      model_preds = model_preds[-4<model_preds$linear_score &  model_preds$linear_score<5,]
    }
    # reg_color
    plot = ggplot(data = model_preds, aes(x=linear_score, color='red')) +
      geom_histogram(data = model_preds, aes(y = ifelse(stat(count) > 0, -stat(count)/(max(model_preds[,SBS])*const), NA)),
                     bins = bins_const, show.legend = FALSE, size=5)+
      geom_point(aes(y = get(SBS), color = "black"), size = 3, show.legend = F)+
      scale_color_manual(values = c(col_SBS1, col_SBS2))+ylim((-1*max(model_preds[,SBS])/8), max(model_preds[,SBS]))+xlim(min(model_preds$linear_score),max(model_preds$linear_score))+
      geom_smooth(aes(y = get(SBS)), method="loess", col = "#7570b3", se = T, linetype = "dashed") +theme_bw()+font("x.text", size = 13) +font("y.text", size = 15) +
      font("xy.title", size = 15) + guides(fill=guide_legend(title="Model")) + theme(strip.text.x = element_text(size=15))+
      labs(x = paste0(tumor_type, ' TP53 GCN Linear Score'), y=paste(SBS))
    tiff(paste0(PATH_figures_regression, "//Regression_", model, "_", tumor_type, "_", gene, "_", SBS, "_with_histogram_3.0.tiff"), units="in", width=6, height=3, res=300)
    print(plot)
    dev.off()


    ### Run boxplot analysis

    model_preds$label = as.factor(model_preds$label)
    levels(model_preds$label) = c("WT", "M")

    model_preds$prediction = as.factor(model_preds$prediction)
    levels(model_preds$prediction) = c("ND", "D")
    # boxplot

    model_preds$pseudo_M_WT = NA
    model_preds$pseudo_M_WT[model_preds$label == "WT" & model_preds$prediction == "ND"] = "tWT"
    model_preds$pseudo_M_WT[model_preds$label == "WT" & model_preds$prediction == "D"] = "pM"
    model_preds$pseudo_M_WT[model_preds$label == "M" & model_preds$prediction == "ND"] = "pWT"
    model_preds$pseudo_M_WT[model_preds$label == "M" & model_preds$prediction == "D"] = "tM"
    model_preds$pseudo_M_WT = factor(model_preds$pseudo_M_WT, c("tWT", "pM", "tM", "pWT"))

    # define side of wilcox.test
    comparisons = list(c("pM", "tWT"), c("tM", "pWT"))
    if(median(model_preds[model_preds$label == "M",SBS]) == 0 & median(model_preds[model_preds$label == "WT",SBS]) == 0) {
      if(mean(model_preds[model_preds$label == "M",SBS]) > mean(model_preds[model_preds$label == "WT",SBS]))  wilc_side = "greater" else wilc_side = "less"
    } else {
      if(median(model_preds[model_preds$label == "M",SBS]) > median(model_preds[model_preds$label == "WT",SBS]))  wilc_side = "greater" else wilc_side = "less"
    }


    # box_label = ggboxplot(model_preds, x= "label", y = SBS,
    #                       color = "label")+labs(x = paste0(gene, " Mutational state\n", tumor_type), y = SBS, fill = "Mutational state")+
    #   stat_compare_means(method = "wilcox.test",size=7, label.y = max(model_preds[,SBS])*0.95)+# method.args = list(alternative = wilc_side))+
    #   geom_jitter(color=brw_set2[3], size=0.4, alpha=0.8)+
    #   scale_color_manual(values = c("#386cb0", "#f0027f"))+font("x.text", size = 20) +font("y.text", size = 20) +
    #   font("xy.title", size = 20) + theme(strip.text.x = element_text(size=20))+ theme(legend.text=element_text(size=20))+ theme(legend.title=element_text(size=20))

    # p_M, p_WT boxplots
     box_pred = ggboxplot(model_preds, x= "pseudo_M_WT", y = SBS,
          color = "pseudo_M_WT")+labs(x = paste0(gene, " predictions\n", tumor_type), y = SBS, fill = "Mutation state")+
    stat_compare_means(comparisons = comparisons, size=7, label.y = max(model_preds[,SBS])*0.91)+#, method.args = list(alternative = wilc_side))+
       geom_jitter(color=brw_set2[3], size=1, alpha=0.8)+scale_color_manual(values = c("#67a9cf", "#1b9e77","#ef8a62","#e6ab02"))
    box_pred = box_pred+guides(color=guide_legend(title=paste0(model, " Predictions")))+font("x.text", size = 22) +font("y.text", size = 25) +
      font("xy.title", size = 25) + theme(strip.text.x = element_text(size=22))+ theme(legend.text=element_text(size=25))+ theme(legend.title=element_text(size=25))

    # model_preds = model_preds[model_preds[,SBS] != 0,]


    tiff(paste0(PATH_figures, "//SBS_boxplots_pseudo_M_WT//p_M_WT//", model, "_", tumor_type, "_", gene, "_", SBS, ".tiff"), units="in", width=10, height=5, res=300)
    print(box_pred)
    dev.off()


    # again with different parameters
    box_label = ggboxplot(model_preds, x= "label", y = SBS,
                          color = "label")+labs(x = paste0(gene, " Mutational state\n", tumor_type), y = SBS, fill = "Mutational state")+
      stat_compare_means(method = "wilcox.test",size=7, label.y = max(model_preds[,SBS])*0.98)+# method.args = list(alternative = wilc_side))+
      geom_jitter(color=brw_set2[3], size=1, alpha=0.8)+
      scale_color_manual(values = c("#67a9cf", "#ef8a62"))+font("x.text", size = 22) +font("y.text", size = 25) +
      font("xy.title", size = 25) + theme(strip.text.x = element_text(size=22))+ theme(legend.text=element_text(size=25))+ theme(legend.title=element_text(size=25))
    # simple pred boxplot
    box_pred = ggboxplot(model_preds, x= "prediction", y = SBS,
                          color = "prediction")+labs(x = paste0("Model predictions\n", tumor_type), y = SBS, fill = "Model predictions")+
      stat_compare_means(method = "wilcox.test", size=7, label.y = max(model_preds[,SBS])*0.98)+#, method.args = list(alternative = wilc_side))+
      geom_jitter(color=brw_set2[3], size=1, alpha=0.8)+
      scale_color_manual(values = c("#67a9cf", "#ef8a62"))+font("x.text", size = 22) +font("y.text", size = 25) +
      font("xy.title", size = 25) + theme(strip.text.x = element_text(size=22))+ theme(legend.text=element_text(size=25))+ theme(legend.title=element_text(size=25))

    tiff(paste0(PATH_figures, "//SBS_boxplots_pseudo_M_WT//pred_label_", model, "_", tumor_type, "_", gene, "_", SBS, ".tiff"), units="in", width=10, height=5, res=300)
    print(box_label + box_pred)
    dev.off()
    
  }
}

# SBS_linear_score_correlations_spearman = cbind(t(as.data.frame(SBS_tumor_type_gene_list)), spearman_rho, spearman_pval)
SBS_linear_score_correlations_spearman = cbind(models, SBS_s, tumor_types, genes, spearman_rho, spearman_pval, wilc_label_pval, wilc_pred_pval)
colnames(SBS_linear_score_correlations_spearman) = c("Model","SBS", "tumor_type", "gene", "spearman_rho", "spearman_pval", "wilc_label_pval", "wilc_pred_pval")
row.names(SBS_linear_score_correlations_spearman) = SBS_linear_score_correlations_spearman[,1]
SBS_linear_score_correlations_spearman = as.data.frame(SBS_linear_score_correlations_spearman)
# generate fdr correction for p values.
pvalues = c(SBS_linear_score_correlations_spearman$spearman_pval,
            SBS_linear_score_correlations_spearman$wilc_label_pval,
            SBS_linear_score_correlations_spearman$wilc_pred_pval)
pvalues_fdr = pvalues_fdr = p.adjust(pvalues, method = "fdr")
SBS_linear_score_correlations_spearman$spearman_pval_fdr = pvalues_fdr[1:16]
SBS_linear_score_correlations_spearman$wilc_label_pval_fdr = pvalues_fdr[17:32]
SBS_linear_score_correlations_spearman$wilc_pred_pval_fdr = pvalues_fdr[33:48]

write.csv(SBS_linear_score_correlations_spearman, paste0(PATH_figures_regression, "//models_SBS_pvalues_wilc_label_pred_spearman_linear_wFDR_3.0.csv"))

