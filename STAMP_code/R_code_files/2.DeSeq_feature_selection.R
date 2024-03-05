
###############################
# create the DeSeq            #
# feature selection process   #
# and spereate train & test   #
# for BRCA                    #
###############################

library(edgeR)
library(dplyr)


PATH = "//media//gil//My Passport//GDCdata//Expression_tables//"
PATH_train = "//home//gil//Desktop//GCN//datasets//BRCA_train_test//"

BRCA_raw_counts_RNA_seq <- read.csv(paste0(PATH, "RNASeq_TCGA_BRCA_raw_counts_for_DeSeq_from_TCGAbiolinks.csv"), row.names = 1)

row.names(BRCA_raw_counts_RNA_seq) <- sapply(row.names(BRCA_raw_counts_RNA_seq), function(x){return(substr(x, 11, 25))})

# 1. load RNA-seq training data (numeric values & labels)
BRCA_train_samples <- read.csv(paste0(PATH_train, "train_Sample_IDs.csv"), header = F)
BRCA_train_label <- read.csv(paste0(PATH_train, "y_train.csv"))
BRCA_train_label <- BRCA_train_label[BRCA_train_samples$V1 %in% row.names(BRCA_raw_counts_RNA_seq),]
BRCA_train_samples <- BRCA_train_samples[BRCA_train_samples$V1 %in% row.names(BRCA_raw_counts_RNA_seq),]

# get only train set's rna seq
BRCA_raw_counts_RNA_seq <- BRCA_raw_counts_RNA_seq[BRCA_train_samples,]


# round all columns
BRCA_raw_counts_RNA_seq <- apply(BRCA_raw_counts_RNA_seq, 2, as.numeric)
# BRCA_raw_counts_RNA_seq <- apply(BRCA_raw_counts_RNA_seq, 2, round)


raw_train_set <- BRCA_raw_counts_RNA_seq

raw_train_label[1:859] <- "ND"
raw_train_label[BRCA_train_label == 1] <- "D"
# use edgeR, limma for differential expression
# input: raw_train_set (only numeric values) and raw_train_label
# output: top positive differentially expressed genes, top negative differentially expressed genes.

# analysis
raw_train_set <- as.data.frame(t(raw_train_set))
d0 <- DGEList(raw_train_set)
d0 <- calcNormFactors(d0)

# some genes are removed
cutoff <- 1
drop <- which(apply(cpm(d0), 1, mean) < cutoff)
d <- d0[-drop,]
# how many genes are left?
# dim(d)

# dimensionality reduction of the samples.
# snames <- colnames(raw_train_set) # Sample names
# plotMDS(d, col = as.numeric(raw_train_label))

# create dif.exp model
mm <- model.matrix(~0+raw_train_label)
y <- voom(d, mm, plot = T)

fit <- lmFit(y, mm)
# head(coef(fit))
contr <- makeContrasts(raw_train_labelD - raw_train_labelND, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
# bayes smoothing
tmp <- eBayes(tmp)

### create table with differentialy expressed genes for D vs. ND. ###
# logFC: log2 fold change of I5.9/I5.6
# AveExpr: Average expression across all samples, in log2 CPM
# t: logFC divided by its standard error
# P.Value: Raw p-value (based on t) from test that logFC differs from 0
top.table <- topTable(tmp, sort.by = "P", n = Inf)
top.table <- top.table[which(top.table$adj.P.Val < 0.05),]
# t is the chosen parameter.
top.table <- top.table[, c("t", "P.Value")]
dif_exp_pos <- top.table[order(top.table$t, decreasing = TRUE),,drop = FALSE]
dif_exp_neg <- top.table[order(top.table$t),,drop = FALSE]

DeSeq_genes_selection <- c(row.names(dif_exp_pos[1:150,]), row.names(dif_exp_neg[1:150,]))
# cut the gene's hugo symbol only
DeSeq_genes_selection <- sapply(DeSeq_genes_selection, function(x){strsplit(x, "[.]")[[1]][1]})
# save
write.csv(DeSeq_genes_selection, row.names = F,
          "//home//gil//Desktop//GCN//datasets//gene_selection_sets//DeSeq_genes_selection_300_set.csv")

