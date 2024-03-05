#############################################
### Create HEATMAP showing                ###
### for each tumor_type and gene/pathway  ###
### the level of learning by test AUC.    ###
#############################################

library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(pheatmap)
library("RColorBrewer")

pathway_path <- "//home//gil//Desktop//GCN//GCN_analysis//6.Cell_pathway_labeling_grid//"
pathway_best_models <- read.csv(paste0(pathway_path, "best_models_GCN_Cell_paper_labels_PATHWAY_FIXED_ADJ.csv"))

# return best auc for tumor type and gene/pathway
get_best_auc <- function(df, tumor_type, gene){
    tumor_type_gene <- df[df$gene == gene & df$cancer.type == tumor_type,]
    max_test_auc <- max(tumor_type_gene$test_auc)
    if(is.infinite(max_test_auc)){return(NA)}
    return(max_test_auc)
}

# order gene list by the maximum num_mutated value available for each gene.
get_gene_list_ordered_by_max_mutated <- function(df){
  order_genes <- data.frame(0)
  for(gene in unique(df$gene)){
    order_genes[gene, "max_num_mutated"] <- max(df$num_mutated[df$gene == gene])
  }
  order_genes <- order_genes[order(order_genes$"max_num_mutated"),]
  return(row.names(order_genes)[-nrow(order_genes)])
}

# order tumor type list by the maximum num_mutated value available for each tumor type
get_tumor_type_list_ordered_by_max_mutated <- function(df){
  order_tumor_types <- data.frame(0)
  for(tumor_type in unique(df$cancer.type)){
    order_tumor_types[tumor_type,"max_num_mutated"] <- max(df$num_mutated[df$cancer.type == tumor_type])
  }
  
  order_tumor_types <- order_tumor_types[order(order_tumor_types$"max_num_mutated"),]
  return(row.names(order_tumor_types)[-nrow(order_tumor_types)])
}

# create table with values for heatmap, using other functions
create_table_for_heatmap <- function(df){
  
  order_genes <- get_gene_list_ordered_by_max_mutated(df)
  order_tumor_types <- get_tumor_type_list_ordered_by_max_mutated(df)
  
  new_df <- data.frame(0)
  for(gene in order_genes){
    for(tumor_type in order_tumor_types){
      new_df[gene,tumor_type] <- get_best_auc(df, tumor_type, gene)
    }
  }
  
  new_df <- new_df[-1,-1]
  
  return(new_df)
}

# create the heatmap figure and save it as "name.tiff" in pathway_path.
create_heatmap_figure <- function(df, name) {
  # use nrow to distinguish gene and ptahway dfs
  if(nrow(df)>150){
    font_row = 30
    font_col = 25
  } else {
    font_row = 20
    font_col = 20
  }
  df_for_heatmap <- create_table_for_heatmap(df)
  
  # create as normalized matrix. NAs get 0.
  matrix <- data.matrix(df_for_heatmap)
  row.names(matrix) <- row.names(df_for_heatmap)
  # matrix_norm[is.na(matrix_norm)] = 0
  
  tiff(paste0(pathway_path, name, ".tiff"), units="in", width=15, height=15, res=300)
  print(pheatmap(matrix, na_col = "grey90",fontsize_row = font_row,fontsize_col = font_col, 
           # scale = "column",      
           cluster_rows = F, cluster_cols = F))
  dev.off()
  
}

############
### MAIN ###
############

# load data
gene_best_models <- read.csv(paste0(pathway_path, "best_models_GCN_Cell_paper_labels_GENE_FIXED_ADJ.csv"))
pathway_best_models <- read.csv(paste0(pathway_path, "best_models_GCN_Cell_paper_labels_PATHWAY_FIXED_ADJ.csv"))
pathway_best_models$gene <- paste0(pathway_best_models$gene, "_pathway")
# pathway_best_models <- pathway_best_models[pathway_best_models$gene %in% c()]
# create heatmaps
create_heatmap_figure(gene_best_models, "Cell_paper_labels_gene_Test_AUC_HEATMAP_3.0")
create_heatmap_figure(pathway_best_models, "Cell_paper_labels_pathway_Test_AUC_HEATMAP_3.0")
