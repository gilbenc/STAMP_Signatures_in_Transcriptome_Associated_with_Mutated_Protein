library(forestmodel)
library(RVAideMemoire)
library(survminer)
library(survival)
library(caret)
library(pROC)
library(dplyr)
library(ggplot2)
library(tidyr)
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
library(mboost)
library(randomForestSRC)
library("RColorBrewer")


### Functions
df_for_survival_pan_cancer_forestplot <- function(pred_table) {
        
        ### load TCGA Clinical data ###
        clinical.path = file.path(paste0(PATH_datasets, "Cbioportal_Pan_Cacner_WmutProfile_clinical_data.tsv"), fsep = .Platform$file.sep)
        clinical_data<- read.csv(clinical.path, sep="\t", fileEncoding="UTF-8-BOM", header=T, stringsAsFactors = F, check.names = F)

        survival_data <- clinical_data[,c('Sample ID', 'Overall Survival (Months)', 'Overall Survival Status',  "TCGA PanCanAtlas Cancer Type Acronym", "Diagnosis Age")]

        # merge with models
        pred_table_survival <- merge(pred_table, survival_data, 'Sample ID')
        
        # survival analysis
        months <- pred_table_survival$`Overall Survival (Months)`
        status <- pred_table_survival$`Overall Survival Status`
        
        label <- pred_table_survival$label
        binary_pred <- pred_table_survival$prediction
        linear_pred <- pred_table_survival$linear_score
        
        tumor_type <- pred_table_survival$`TCGA PanCanAtlas Cancer Type Acronym`
        diagnosis_age <- pred_table_survival$`Diagnosis Age`
        # categorical_pred <- # to be modified
        
        status <- as.integer(status == "1:DECEASED")
        dataframe <- data.frame(months, status, binary_pred, linear_pred, label, tumor_type, diagnosis_age)
        return(dataframe)
}

### Load pan cancer model (using 'Survival_analysis_and_curves.R')
load("~/Desktop/GCN/Main_code_scripts/R_scripts/R_workspaces/survival_model_pan_cancer_for_forest_plot.RData")


df_merged_survival = df_for_survival_pan_cancer_forestplot(pred_table)
colnames(df_merged_survival)

# remove 3 tumor types with too high variance (obstruct the others. contain very small number of samples (PCPG, TGCT: n =1, DLBC: n = 3)
# mltvar_data_survival <- mltvar_data_survival[!mltvar_data_survival$tumor_type %in% c("PCPG", "TGCT", "DLBC"),]
df_merged_survival$tumor_type <- factor(df_merged_survival$tumor_type, levels = c("LUAD", unique(df_merged_survival$tumor_type[-which(df_merged_survival$tumor_type== "LUAD")])))
df_merged_survival$binary_pred[df_merged_survival$binary_pred == 0] <- "ND"
df_merged_survival$binary_pred[df_merged_survival$binary_pred == 1] <- "D"
df_merged_survival$binary_pred <- factor(df_merged_survival$binary_pred)
levels(df_merged_survival$binary_pred) = c("ND", "D")
# fit <- coxph(Surv(mltvar_data_survival$OS_MONTHS, mltvar_data_survival$OS_STATUS)~mltvar_data_survival$TP53_PROF+mltvar_data_survival$tumor_type)
fit <-coxph(Surv(months, status)~binary_pred+tumor_type, data = df_merged_survival)

### until here: load a survival model with the variables you wish to examine.
load("~//Desktop//GCN//Main_code_scripts//R_scripts//R_workspaces//RANK_survival_model_pan_cancer_for_forest_plot_09_22.RData")
dataframe$quantiles_label <- quantiles_label
fit_pan_cancer_multivariable <- coxph(Surv(months, status) ~ factor(quantiles_label)+tumor_type_col+age+Mutation_count, data = dataframe)
fit_pan_cancer_multivariable_quartiles_as_factor <- coxph(Surv(months, status) ~ factor(quantiles_label)+tumor_type_col+age+Mutation_count, data = dataframe)

forest <- forest_model(fit_pan_cancer_multivariable)
forest_q_as_factor <- forest_model(fit_pan_cancer_multivariable_quartiles_as_factor)


tiff("//home//gil//Desktop//GCN//GCN_analysis//Figures_for_paper//Figure_4_Survival//survival_forestplot//RANKING_by_Quartiles_pan_cancer_forestplot_multivariate_tumor_type_age_mut_burden.tiff", units="in", width=8, height=3, res=200)
print(plot(forest))
dev.off()

tiff("//home//gil//Desktop//GCN//GCN_analysis//Figures_for_paper//Figure_4_Survival//survival_forestplot//RANKING_by_Quartiles_AS_FACTOR_pan_cancer_forestplot_multivariate_tumor_type_age_mut_burden.tiff", units="in", width=8, height=3, res=200)
print(plot(forest_q_as_factor))
dev.off()


plot(forest)

forest_df <- data.frame(forest$data$level, forest$data$estimate, forest$data$conf.low, forest$data$conf.high, forest$data$n, stringsAsFactors = F)

forest_df[c(1,5),c(3,4)] <- 0
forest_df[36,1] <- "Age"
forest_df[37,1] <- "Mutation Burden"
forest_df[1:4,"setting"] <- "STAMP Quartiles"
forest_df[c(5:35), "setting"] <- "Tumor Type"
forest_df[c(36), "setting"] <- "Age"
forest_df[c(37), "setting"] <- "Mut"

# forest_df[36,c("forest.data.level","setting")] <- "Diagnosis Age"
colnames(forest_df) <- c("Reference", "effect_size", "CI_low", "CI_high", "n", "setting")
forest_df[,1] <- paste0(forest_df$Reference, " (n = ", forest_df$n, ")")
forest_df[c(1,5),1] <- paste0(forest_df[c(1,5),1], " Ref")
forest_df[,1] <- factor(forest_df[,1], levels = rev(forest_df[,1]))
forest_df$setting <- factor(forest_df$setting, levels = unique(forest_df$setting))
# add column for distinct p-val tumor types
forest_df$distinct[forest$data$p.value < 0.05 & forest_df$setting != "Pred"] <- 1
forest_df$pval <- forest$data$p.value

# Make a plot called 'p', and map citation data to y-axis, effect sizes to x-axis
# specify the min and max of the CIs, and give different shapes based on levels of tester
p=ggplot(forest_df, aes(y=Reference, x=exp(effect_size), xmin=exp(CI_low), xmax=exp(CI_high)))+
        #Add data points and color them black
        geom_point(data=subset(forest_df, setting!= "Pred"),color = 'black', size = 3)+
        #add the CI error bars
        geom_errorbarh(height=.1, size = 1)+
        #Specify the limits of the x-axis and relabel it to something more meaningful
        scale_x_continuous(limits=c(0.01,7), trans = "log", name='Hazard Ratio', breaks = c(0.1,0.2, 0.5, 1,2,3, 5), labels=c("0.1","0.2","0.5","1", "2", "3", "5"))+
        #Give y-axis a meaningful label
        ylab('')+
        #Add a vertical dashed line indicating an effect size of zero, for reference
        geom_vline(xintercept=1, color='black', linetype='dashed') + 
        theme_bw() +
        facet_grid(setting~., scales='free', space='free') +
        theme(axis.text.x =element_text(size=10, face = "bold")) +
        theme(axis.text.y =element_text(size=10, face = "bold", family = "sans")) +
        #Add 'special' points for the summary estimates, by making them diamond shaped
        # geom_point(data=subset(forest_df, Reference=='D (n = 5685)'), color='#80B1D3', shape=18, size=4)+
        # geom_errorbarh(data = subset(forest_df,Reference=='D (n = 5685)'), color='#80B1D3', height=.1, size = 1)+
        # geom_point(data=subset(forest_df, Reference=='ND (n = 3158) Ref'), color='#FFED6F', shape=18, size=4)
        geom_point(data=subset(forest_df, distinct==1), color='#E18F15', size=3)+
        geom_errorbarh(data = subset(forest_df,distinct==1), color='#E18F15', height=.1, size = 1)


tiff("//home//gil//Desktop//GCN//GCN_analysis//Figures_for_paper//Figure_4_Survival//survival_forestplot//RANKING_by_Quartiles_AS_FACTOR_pan_cancer_forestplot_multivariate_tumor_type_age_mut_burden_GGPLOT.tiff", units="in", width=8, height=10, res=200)
print(plot(p))
dev.off()


print(forest_df$pval[2])

tiff("//home//gil//Desktop//GCN//GCN_analysis//Figures_for_paper//Figure_4_Survival//survival_forestplot//forestplot_multivariate_tumor_type.tiff", units="in", width=5, height=6, res=200)
print(p)
dev.off()
