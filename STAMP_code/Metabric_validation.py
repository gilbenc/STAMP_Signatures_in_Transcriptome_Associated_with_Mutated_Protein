import pandas as pd
import torch
import os, fnmatch
import numpy as np
from data.datasets import TCGADataset
from matplotlib import pyplot as plt
import scipy
import math
import pickle
import Parse_data as Parse
from data.gene_graphs import GeneManiaGraph, RegNetGraph, HumanNetV2Graph, \
    FunCoupGraph
from data.datasets import TCGADataset

### Functions

# returns model of model_type, using data_type, for cancer_type & gene combination
def load_model(data_type, model_type, cancer_type, gene):
    PATH_model = "/home/gil/Desktop/GCN/GCN_analysis/2.models/" + data_type + "/" + model_type + "/"
    if model_type == "RF":
        pattern = "*" + '_' + cancer_type + "_" + gene + "_deseq_" + "*"
    else:
        pattern = "*" + '_' + cancer_type + "_" + gene + "*"

    model_path = find(pattern, PATH_model)
    if len(model_path) > 1:
        if model_type == "ELR":
            pattern = "*" + '_' + cancer_type + "_" + gene + "_mutsigdb_" + "*"
            model_path = find(pattern, PATH_model)
            if len(model_path) > 1:
                print("error, too many files fit description.")
                exit()
        else:
            print("error, too many files fit description.")
            exit()

    if model_type == "GCN":
        return torch.load(model_path[0])
    else:
        # filename = PATH + "GCN_analysis/2.models/TP53_many_tumor_types/" + model_type + "/best_model_" + model_type + "_" + cancer_type + "_" + gene + "_" + gene_selection + "_AUC_" + str(
        #     best_auc) + ".sav"

        return pickle.load(open(model_path[0], 'rb'))


# subfunction for load_model
# return 'result' containing file name (should be only one) in 'path' with 'pattern' in them
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

# Load GDSC rnaseq data and parse:
# a. keep only 1st degree
# b. change column names to cell line names
# c. get rid of unnecessary columns
# d. transpose
# e. compute zscore(log(rsem values))
# f. add first degree missing columns as 0 values
def parse_METABRIC_rnaseq(first_degree_neighbors, cancer_type):
    data_PATH = "/home/gil/Desktop/GCN/datasets/Data_downloaded_directly_from_source/"

    # read METABRIC rnaseq
    METABRIC_df = pd.read_csv(data_PATH + 'METABRIC_gene_expression_ILMN.csv')
    # read GDSC cell line-target-drugs data data
    # print("now what?")

    # COSMIC_to_Cell_line_name = GDSC_cell_target_drug.loc[:, ['COSMIC_ID', 'CELL_LINE_NAME']].drop_duplicates()
    # COSMIC_to_Cell_line_name.index = COSMIC_to_Cell_line_name['COSMIC_ID']

    ### Do I have hugo symbols???
    # get hugo values in intersect of GDSC and 1st deg neighbors
    in_first_deg_and_METABRIC = np.unique([x for x in METABRIC_df['GENE_SYMBOLS'] if x in first_degree_neighbors])
    # use hugo as index so I can call loc for 1st deg neighbors.
    METABRIC_df.index = METABRIC_df['GENE_SYMBOLS']
    # keep in METABRIC only intersect genes, by the order of the intersect vector
    METABRIC_df = METABRIC_df.loc[in_first_deg_and_METABRIC]
    # remove duplicate gene symbols: simply take first one (assuming values are similar for genes)
    METABRIC_df = METABRIC_df[~METABRIC_df.index.duplicated(keep='first')]



    ### I will need sample IDs for this one
    # replace GDSC_df's columns with cell_line_names, transpose, remove duplicates
    # GDSC_rna_cosmic_values = list(set([int(str(col).split(".")[1]) for col in GDSC_df.columns[2:]]))
    # GDSC_rna_cosmic_values = [val for val in GDSC_rna_cosmic_values if val in COSMIC_to_Cell_line_name.index]
    # COSMIC_to_Cell_line_name = COSMIC_to_Cell_line_name.loc[GDSC_rna_cosmic_values]
    # # Transpose, get rid of first 2 columns
    # GDSC_df = GDSC_df.iloc[:,2:].transpose()
    # GDSC_df.index = [int(str(ind).split(".")[1]) for ind in GDSC_df.index]
    # GDSC_df = GDSC_df.loc[GDSC_rna_cosmic_values].drop_duplicates()
    # GDSC_df.index = COSMIC_to_Cell_line_name.loc[GDSC_df.index, 'CELL_LINE_NAME']

    # REMOVE irrelveant columns first/last and then transpose.
    metabric_sample_IDs = METABRIC_df.columns
    METABRIC_df = METABRIC_df.iloc[:,1:2115].transpose()


    # Parse METABRIC into z scores similar to TCGA - zscore(log(rna_seq values))
    METABRIC_df = np.log(METABRIC_df + (1e-10))
    for col in METABRIC_df.columns:
        METABRIC_df[col] = scipy.stats.zscore(METABRIC_df[col])


    # add and reorder METABRIC columns to fit TCGA
    add_columns = [col for col in first_degree_neighbors if col not in METABRIC_df.columns]
    print("need to add "+ str(add_columns) + " to METABRIC gene exp. to fit TCGA.")
    for col in add_columns:
        METABRIC_df[col] = 0
    METABRIC_df = METABRIC_df.reindex(columns=first_degree_neighbors)

    return METABRIC_df

if __name__ == '__main__':
    # paths
    PATH_GCN = "/home/gil/Desktop/GCN/"
    data_PATH = "/home/gil/Desktop/GCN/datasets/Data_downloaded_directly_from_source/"
    results_PATH = "/home/gil/Desktop/GCN/GCN_analysis/11.Metabric_validation/"

    # define model/gene/cancer_type/drug to investigate
    # model_type = "GCN"
    delay_seconds = 5
    cancer_type = "BRCA"
    use_pan_cancer = False
    gene = "TP53" # or ERBB2 !!
    for model_type in ["ELR"]: # ['GCN', 'ELR', 'RF']

        # Also define first deg neighbors (GCN & ELR. for RF its actually deseq genes selction)
        if model_type == "GCN":


            # get model
            if not use_pan_cancer:
                data_type = "Pathways"
                gene_model = load_model(data_type, model_type, cancer_type, gene)
            else:
                data_type = "TP53_many_tumor_types"
                gene_model = load_model(data_type, model_type, "pan_cancer", gene)

            gene_first_degree_neighbors = gene_model.X.columns.drop_duplicates()

        else:
            # Read in data: TCGA
            dataset = TCGADataset()
            # Remove duplicate genes!
            dataset.df = dataset.df.loc[:, ~dataset.df.columns.duplicated()]
            # save a copy of dataset, since it will be modified for cancer_type specific samples
            dataset_df_copy = dataset.df

        if model_type == "RF":
            data_type = "TP53_many_tumor_types"  # "Pathways"
            gene_set = Parse.get_genes_selection_set("deseq")
            gene_first_degree_neighbors = [n for n in gene_set if n in dataset_df_copy.columns.values]

            # get model
            gene_model = load_model(data_type, model_type, cancer_type, gene)
        if model_type == "ELR":
            data_type = "TP53_many_tumor_types"
            gene_set = Parse.get_genes_selection_set("mutsigdb")
            gene_first_degree_neighbors = [n for n in gene_set if n in dataset_df_copy.columns.values]

            # get model
            if not use_pan_cancer:
                gene_model = load_model(data_type, model_type, cancer_type, gene)
            else:
                gene_model = load_model(data_type, model_type, "pan_cancer", gene)
            # gene_graph_funcoup = FunCoupGraph()
            # neighbors_funcoup = list(gene_graph_funcoup.first_degree(gene)[0])
            # gene_first_degree_neighbors = [n for n in neighbors_funcoup if n in dataset_df_copy.columns.values]




        # 1st degree neighbors


        # parsre METABRIC for model predictions:
        METABRIC_first_deg_z_scored = parse_METABRIC_rnaseq(gene_first_degree_neighbors, cancer_type)
        # Impute if necessary
        if METABRIC_first_deg_z_scored.isnull().values.any():
            METABRIC_first_deg_z_scored.fillna(0)

        # get predictions for Metabric
        gene_model_predict = gene_model.predict(METABRIC_first_deg_z_scored)

        if model_type == "GCN":
            Metabric_predict_bool = np.argmax(gene_model_predict, axis=1)
            # calculate linear score, avoid infinity values
            Metabric_predict_linear = scipy.special.logit(gene_model_predict[:, 1])
            max_value = max(Metabric_predict_linear[Metabric_predict_linear != math.inf])
            min_value = min(Metabric_predict_linear[Metabric_predict_linear != -math.inf])
            Metabric_predict_linear[Metabric_predict_linear == math.inf] = max_value + 1
            Metabric_predict_linear[Metabric_predict_linear == -math.inf] = min_value - 1

            # create Metabric dataframe with GCN preds
            Metabric_predictions = pd.DataFrame(METABRIC_first_deg_z_scored.index)
            Metabric_predictions.index = Metabric_predictions.iloc[:, 0]
            Metabric_predictions = Metabric_predictions.rename(columns={0: "Metabric_Sample_ID"})

            Metabric_predictions[model_type + "_pred"] = Metabric_predict_bool.numpy()
            Metabric_predictions[model_type + "_linear"] = Metabric_predict_linear.numpy()
            Metabric_predictions = Metabric_predictions.iloc[:, 1:]

        else:
            if model_type == "RF":
                Metabric_predict_bool = gene_model_predict
                Metabric_predict_linear = gene_model.predict_proba(METABRIC_first_deg_z_scored)[:, 0]
            else:
                Metabric_predict_bool = (gene_model_predict > 0.5).astype(int)
                # calculate linear score, avoid infinity values
                Metabric_predict_linear = gene_model_predict
            # create Metabric dataframe with GCN preds
            Metabric_predictions = pd.DataFrame(METABRIC_first_deg_z_scored.index)
            Metabric_predictions.index = Metabric_predictions.iloc[:, 0]
            Metabric_predictions = Metabric_predictions.rename(columns={0: "Metabric_Sample_ID"})

            Metabric_predictions[model_type + "_pred"] = Metabric_predict_bool
            Metabric_predictions[model_type + "_linear"] = Metabric_predict_linear
            Metabric_predictions = Metabric_predictions.iloc[:, 1:]

        if not use_pan_cancer:
            Metabric_predictions.to_csv(
                results_PATH + "METABRIC_zscored_pan_" + gene + "_" + cancer_type + "_" + model_type + "_predictions.csv"
            )
        else:
            Metabric_predictions.to_csv(
                results_PATH + "METABRIC_zscored_pan_" + gene + "_PAN_CANCER_" + model_type + "_predictions.csv"
            )
        del gene_model
        if model_type == "GCN":
            torch.cuda.empty_cache()


