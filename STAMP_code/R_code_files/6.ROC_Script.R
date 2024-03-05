### ROC Script for Natan ###

library(ggplot2)
library(ggpubr)
library("RColorBrewer")
library(ROCR)
library(tidyr)
library(pROC)



### Load Data


## A. create performance for ROC curve calculations
# label should be metabric's label for mutated/non-mutated TP53, pred_score is the model's prediction score
ELR_pred <- prediction(ELR_pred_score, label)
ELR_perf <- performance(pred, "tpr", "fpr")

GCN_pred <- prediction(GCN_pred_score, label)
GCN_perf <- performance(pred, "tpr", "fpr")

RF_pred <- prediction(RF_pred_score, label)
RF_perf <- performance(pred, "tpr", "fpr")


## B. create figure
tiff(paste0(PATH_plots, "METABRIC_ROC_3_Models_qualitative_3.0.tiff"),
     units="in", width=8, height=8, res=300)
par(cex.axis=2, cex.lab = 2, mar=c(5,6,4,1))
print(plot(get(paste0("perf_", tumor_type, "_GCN")), 
           col = brw_dark2[1], lwd = 7)) 
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