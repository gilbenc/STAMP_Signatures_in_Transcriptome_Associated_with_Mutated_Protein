import pandas as pd
import torch
import os, fnmatch
import numpy as np
import scipy
import math
import pickle
import warnings
from datetime import datetime

### Functions

# returns model of model_type, using data_type, for cancer_type & gene combination
def load_model(data_type, model_type, cancer_type, gene):
    PATH_model = "/models/" + data_type + "/" + model_type + "/"
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

# Load input rnaseq data and parse:
# a. keep only 1st degree
# b. change column names to cell line names
# c. get rid of unnecessary columns
# d. transpose
# e. compute zscore(log(rsem values))
# f. add first degree missing columns as 0 values
def parse_input_rnaseq(first_degree_neighbors):

    # read input RNA data
    input_gene_expression_df = pd.read_csv('gene_expression_input.csv')

    ### Do I have hugo symbols?
    # get hugo values in intersect of input data and 1st deg neighbors
    in_first_deg_and_input = np.unique([x for x in input_gene_expression_df['GENE_SYMBOLS'] if x in first_degree_neighbors])
    # use hugo as index so I can call loc for 1st deg neighbors.
    input_gene_expression_df.index = input_gene_expression_df['GENE_SYMBOLS']
    # keep in input only intersect genes, by the order of the intersect vector
    input_gene_expression_df = input_gene_expression_df.loc[in_first_deg_and_input]
    # remove duplicate gene symbols: simply take first one (assuming values are similar for genes)
    input_gene_expression_df = input_gene_expression_df[~input_gene_expression_df.index.duplicated(keep='first')]



    # REMOVE irrelveant columns first/last and then transpose.
    input_sample_IDs = input_gene_expression_df.columns
    input_gene_expression_df = input_gene_expression_df.iloc[:,1:2115].transpose()


    # Parse input into z scores similar to TCGA - zscore(log(rna_seq values))
    input_gene_expression_df = np.log(input_gene_expression_df + (1e-10))
    for col in input_gene_expression_df.columns:
        input_gene_expression_df[col] = scipy.stats.zscore(input_gene_expression_df[col])


    # add and reorder input columns to fit TCGA
    add_columns = [col for col in first_degree_neighbors if col not in input_gene_expression_df.columns]
    print("\n"+str(len(add_columns)) + " genes are missing from input to fit TCGA. columns are added with 0 values. This may affect predictions.")
    print("genes added with 0 values are: \n\n" + str(add_columns))
    for col in add_columns:
        input_gene_expression_df[col] = 0
    input_gene_expression_df = input_gene_expression_df.reindex(columns=first_degree_neighbors)

    return input_gene_expression_df

if __name__ == '__main__':

    # Enter file names here for data files you want to run STAMP on
    data_list = []
    # Enter model file names here for the models you wish to apply
    ### IN THE SAME ORDER AS THE DATA FILES THEY SHOULD BE APPLIED TO ###
    models_list = []
    
    
    mrna_data_path = "/mrna_data/"
    STAMP_models_path = "/STAMP_models/"

    for data in data_list:
        
    # Do you want to use one model file (Simple) or multiple files for different tumor types/genes/pathways (Advanced)
    while True:
        input = int(input("Use one model file (Simple) or multiple files for different tumor types/genes/pathways (Advanced)? (type: 1 for simple; 2 for advanced)\n"))
        if not(input in [1, 2]):
            print("input is invalid (type: 1 for simple; 2 for advanced)")
            continue
        else:
            if input == 1:
                input_type = "simple"
            if input == 2:
                input_type = "advanced"
            ...
            break

    if input_type == "simple":
        print("\nMake sure the model file is in the same directory as this script, named 'model.pt'")
        print("Regarding input gene expression data file, ensure that:")
        print("     1. It is in the same directory as this script")
        print("     2. It is a csv file named 'gene_expression_input.csv'")
        print("     3. It has gene names (Hugo symbols) in a column named 'GENE_SYMBOLS' and samples as columns")
        print("        (The script will transpose this to generate predictions. first row of csv should be sample ids/names)")
    else:
        print("\nMake sure the model files are in the /model/data_type/model_type/cancer_type/ directories \n"
              "as ordered by relevant input values, and named with the relevant gene/pathway, ending with '.pt'")
        print("Regarding input gene expression data file/s, ensure that:")
        print("     1. It is in the same directory as this script")
        print("     2.  It is a csv file named 'gene_expression_input.csv'")
        print("     3. It has gene names (Hugo symbols) in a column named 'GENE_SYMBOLS' and samples as columns")
        print("        (The script will transpose this to generate predictions. first row of csv should be sample ids/names)")


    # Check for results directory, if doesn't exist create it
    results_PATH = ("results/")
    CHECK_FOLDER = os.path.isdir(results_PATH)
    if not CHECK_FOLDER:
        os.makedirs(results_PATH)
        print("\ncreated folder : ", results_PATH)


    # meta variables
    use_pan_cancer = False
    data_type = "Pathways"
    model_type = "GCN"

    if input_type == "simple":
        gene_model = torch.load("model.pt")
    else:
        # todo: remove these after applying input from user
        cancer_type = "BRCA" # remove this once input is corrected
        gene = "TP53" # remove this once input is corrected

        # todo: query for tumor type, gene/pathway.
        #  Write script to allow running models based on tumor type and gene if you want to run multiple models on after the other (instead of user input)
        #  raise error and require new input if there's no model for requested query

        # gene/pathway selection by input:
        # todo Modify
        # while True:
            #gene_input = input("Which gene/pathway to use? (type: 1 for ; 2 for ; 3 for ...)")
            #if not(gene_input in [1:?]):
                #print("input is invalid (type: 1 for ; 2 for ;...)")
            #     continue
            # else:
                # if gene_input == 1:
                #     gene = ""
                # if gene_input == 2:
                #     gene = ""
                # ...
                # break

        # cancer_type selection by input:
        # todo Modify
        # while True:
            #cancer_type_input = input("Which cancer_type is your data? (type: 1 for ; 2 for ; 3 for ...)")
            #if not(cancer_type_input in [1:?]):
                #print("input is invalid (type: 1 for ; 2 for ;...)")
            #     continue
            # else:
                # if cancer_type_input == 1:
                #     cancer_type = ""
                # if cancer_type_input == 2:
                #     cancer_type = ""
                # ...
                # break



        # Get model
        gene_model = load_model(data_type, model_type, cancer_type, gene)

    # After model was loaded:
    # Define first deg neighbors
    gene_first_degree_neighbors = gene_model.X.columns.drop_duplicates()



    # parsre input for model predictions:
    input_first_deg_z_scored = parse_input_rnaseq(gene_first_degree_neighbors)
    # Impute if necessary
    if input_first_deg_z_scored.isnull().values.any():
        input_first_deg_z_scored.fillna(0)

    # get predictions for input
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gene_model_predict = gene_model.predict(input_first_deg_z_scored)


    input_predict_bool = np.argmax(gene_model_predict, axis=1)
    # calculate linear score, avoid infinity values
    input_predict_linear = scipy.special.logit(gene_model_predict[:, 1])
    max_value = max(input_predict_linear[input_predict_linear != math.inf])
    min_value = min(input_predict_linear[input_predict_linear != -math.inf])
    input_predict_linear[input_predict_linear == math.inf] = max_value + 1
    input_predict_linear[input_predict_linear == -math.inf] = min_value - 1

    # create input dataframe with GCN preds
    input_predictions = pd.DataFrame(input_first_deg_z_scored.index)
    input_predictions.index = input_predictions.iloc[:, 0]
    input_predictions = input_predictions.rename(columns={0: "input_Sample_ID"})

    input_predictions["STAMP_pred"] = input_predict_bool.numpy()
    input_predictions["STAMP_linear"] = input_predict_linear.numpy()
    input_predictions = input_predictions.iloc[:, 1:]

    now = datetime.now()

    input_predictions.to_csv(
        results_PATH + "input_" + model_type + "_predictions_"+now.strftime("%Y_%m_%d_%H_%M")+".csv"
    )

    print("\nSTAMP predictions saved to: ")
    print("results/input_" + model_type + "_predictions_"+now.strftime("%Y_%m_%d_%H_%M_%S")+".csv")

    del gene_model
    if model_type == "GCN":
        torch.cuda.empty_cache()


