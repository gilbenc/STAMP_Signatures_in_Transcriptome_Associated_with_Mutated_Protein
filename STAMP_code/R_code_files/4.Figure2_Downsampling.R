
library(caret)
library(pROC)
library(dplyr)
library(glmnet)
library(tibble)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(glue)
#library(ModelMetrics)
library(OpenMPController) # for Kaggle backend
library(readr)
library(vtreat)
library(xgboost)
library(ComplexHeatmap)
library(reshape)
library(ggfortify)
library(ramify)
library(e1071)
library(Matrix)
library(DiagrammeR)
library(reshape2)
library(igraph)
library(plotly)



### Preapare table with AUC, F1 ###

PATH_reduction_plot <-"//home//gil//Desktop//GCN//GCN_analysis//3.Models_comparison_grid//"
sample_reduction_figure <- read.csv(paste0(PATH_reduction_plot, 
                                           "best_models_BRCA_fixed_TEST_sample_reduction_LR_ELR_RF_MLP_GCN.csv"))
# Downsampling values
sample_reduction_figure$mut_percent <- 0.32 - sample_reduction_figure$percent_reduce
sample_reduction_figure$n_samples <- 600 - sample_reduction_figure$num_reduce

# F1
test_prec <- sample_reduction_figure$test_prec
test_recall <- sample_reduction_figure$test_recall
sample_reduction_figure$test_F1 <- 2*((test_prec*test_recall)/(test_prec+test_recall))
# Wherever prec&recall = 0, it seems the model simply predicted all 0, therefore, F1 gets 0 as well.
# (test accuracy consistently = % of ND samples in test set)
sample_reduction_figure$test_F1[is.nan(sample_reduction_figure$test_F1)] <- 0

sample_reduction_figure <- sample_reduction_figure[sample_reduction_figure$model != "LR",c("model", "gene_selection", "cancer.type", "n_samples", 
                                                      "mut_percent", "test_auc", "test_F1")]


### Selected approach: Heatmap ###
GCN_AUC <-sample_reduction_figure[sample_reduction_figure$model == "GCN", c("n_samples", "mut_percent", "test_auc")]
GCN_F1 <- sample_reduction_figure[sample_reduction_figure$model == "GCN", c("n_samples", "mut_percent", "test_F1")]
ELR_AUC <- sample_reduction_figure[sample_reduction_figure$model == "ELR", c("n_samples", "mut_percent", "test_auc")]
ELR_F1 <- sample_reduction_figure[sample_reduction_figure$model == "ELR", c("n_samples", "mut_percent", "test_F1")]

GCN_matrix <- pivot_wider(GCN_F1, names_from = n_samples, values_from = test_F1)
GCN_matrix <- as.matrix(GCN_matrix)
row.names(GCN_matrix) <- GCN_matrix[,1]
GCN_matrix <- GCN_matrix[,-1]

ELR_matrix <- pivot_wider(ELR_F1, names_from = n_samples, values_from = test_F1)
ELR_matrix <- as.matrix(ELR_matrix)
row.names(ELR_matrix) <- ELR_matrix[,1]
ELR_matrix <- ELR_matrix[,-1]

PATH_figures <- "//home//gil//Desktop//GCN//GCN_analysis//Figures_for_paper//Figure_3_Downsampling_and_all_tumor_types//"
tiff(paste0(PATH_figures, "ELR_AUC_heatmap.tiff"), units="in", width=8, height=8, res=300)
ggplot(GCN_F1,aes(factor(n_samples, levels = seq(100, 600, 100)),factor(mut_percent, levels = seq(0.05, 0.3, 0.05)))) + 
        geom_tile(aes(fill=test_F1),color = "white") +
        guides(fill=guide_colorbar("AUC Score")) +
        scale_fill_gradientn(colors = c("#2c7bb6", "#ffffbf", "#d7191c"), guide="colorbar",
                             breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2)) +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + 
        theme(axis.title.x = element_text(angle = 0, vjust = -0.3)) + 
        theme(axis.title.y = element_text(vjust = 2)) + 
        theme(legend.title = element_text(vjust = 2,  size = 20)) +
        theme(legend.text = element_text(size = 20)) +
        font("xy.text", size = 20) +  font("xy.title", size = 20) +  
        labs(x = "Number Of Samples", y = "% Mutated Samples") + expand_limits(fill = c(0, 0.9))
dev.off()
# print(p)


### Option 2: similar to tumor types plot (categorical scatter plot+regression line)


# select the variable: score
score = "F1"
# score = "AUC"
downsample <- "n_samples"
# downsample <- "mut_percent"
brw_dark2 <- brewer.pal(n = 3, name = "Dark2")

tiff(paste0("//home//gil//Desktop//GCN//GCN_analysis//Figures_for_paper//Figure_3_Downsampling_and_all_tumor_types///Down_sampling_models_comparison_by_", downsample, "_", 
           score, ".tiff"), units="in", width=24, height=8, res=300)
p <- ggplot(data = sample_reduction_figure,
            aes(x = get(downsample), y = get(paste0("test_", if(score == "AUC") tolower(score) else score)),  colour = model, shape = model)) +
        geom_point(size = 10) + labs(if(downsample == "mut_percent") x = "Mutated Percentage" else x = "Number Of Samples", 
                                     y = paste0("Test ", score)) + theme_bw()+
        stat_smooth(method = "loess", se = F, size = 5) + scale_color_manual(values = brw_dark2)+
        theme(plot.subtitle = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        font("x.text", size = 30) +font("y.text", size = 30) +
        font("legend.title", size = 30) + font("legend.text", size = 30) + 
        font("xy.title", size = 30) #+ labs(color=l, shape=l)#+ guides(col=guide_legend(title="Models")) #+
# scale_x_discrete(labels=unique(test_AUC_models_tumor_types$tumor_type))
print(p)
dev.off()


### Option 3: 3D ###
PATH_reduction_plot <-"//home//gil//Desktop//GCN//GCN_analysis//3.Models_comparison_grid//"
sample_reduction_figure_3D <- read.csv(paste0(PATH_reduction_plot, "best_models_BRCA_fixed_TEST_sample_reduction_LR_ELR_RF_MLP_GCN.csv"))

# calculate F1
test_prec <- sample_reduction_figure_3D$test_prec
test_recall <- sample_reduction_figure_3D$test_recall
sample_reduction_figure_3D$test_F1 <- 2*((test_prec*test_recall)/(test_prec+test_recall))

# calculate new_score
test_auc <- sample_reduction_figure_3D$test_auc
test_percentage <- 0.32 - sample_reduction_figure_3D$percent_reduce
sample_reduction_figure_3D$test_new_score <- (test_auc - (1-test_percentage))/test_percentage
sample_reduction_figure_3D$test_new_score[sample_reduction_figure_3D$test_new_score < 0] <- 0



sample_reduction_figure_3D$test_F1[is.na(sample_reduction_figure_3D$test_F1)] <- 0.3
GCN_df <- sample_reduction_figure_3D[sample_reduction_figure_3D$model == "GCN",]
ELR_df <- sample_reduction_figure_3D[sample_reduction_figure_3D$model == "ELR",]
RF_df <- sample_reduction_figure_3D[sample_reduction_figure_3D$model == "RF",]
LR_df <- sample_reduction_figure_3D[sample_reduction_figure_3D$model == "LR",]



create_3D_surface_plot <- function(model_type, metric) {
        
        model_df <- get(paste0(model_type, "_df"))
        df <- data.frame(x = as.numeric(600 - model_df$num_reduce), 
                         y = as.numeric(0.32 - model_df$percent_reduce),
                         z = as.numeric(model_df[,metric]))
        g=graph.data.frame(df,directed=FALSE)
        matrix <-dcast(data = df, formula = x~y, fun.aggregate =  mean)
        row.names(matrix) <- matrix[,1]
        matrix <- matrix[,-1]
        
        # set bar labels
        f <- list(family = "Arial",
                  size = 18,color = "#7f7f7f")
        x <- list(title = "TP53 Mutated Percentage",titlefont = f)
        y <- list(title = "Number of Samples",titlefont = f)
        test_AUC <- data.matrix(matrix)
        plot_ly(x =  ~colnames(matrix), y = ~row.names(matrix),z = ~test_AUC)%>%
                add_surface(colorscale = list(c(0, 0.6, 1), c("#80B1D3", "#8DD3C7", "#FB8072"))) %>% layout(
                                title = model_type,
                                
                                scene = list(
                                        xaxis = list(title = "TP53 Mutated Percentage", showspikes=FALSE),
                                        yaxis = list(title = "Number of Samples", showspikes=FALSE),
                                        zaxis = list(title = paste0(metric), showspikes=FALSE)
                                )
                        )
}

create_3D_surface_plot("GCN", "test_F1")
create_3D_surface_plot("ELR", "test_F1")
create_3D_surface_plot("RF", "test_F1")

create_3D_surface_plot("GCN", "test_new_score")
create_3D_surface_plot("ELR", "test_new_score")
create_3D_surface_plot("RF", "test_new_score")

create_3D_surface_plot("GCN", "test_acc")
create_3D_surface_plot("ELR", "test_acc")
create_3D_surface_plot("RF", "test_acc")

