To run STAMP on new gene expression data with selected model (tumor type and gene/pathway):
 
    1. Get files from 'use_STAMP_on_new_data' directory in this git.
 
    2. Install STAMP’s environment using anaconda.
        a. If not installed - install anaconda
        b. Env file: Install STAMP_evn.yml as a conda environment 
	   (type in terminal: ‘conda env create -f STAMP_env.yml’ ). 
	c. Activate the environment (type in terminal: ‘conda activate STAMP_env’)

    3. Keep these files in the same directory (navigate command line to this directory):
        a. Script file: named “Run_STAMP_by_user_input.py” .
        b. Model file (Download/receive from STAMP creators); 
	   Name the model file “model.pt” (or "model.sav" if its ELR/RF).
        c. Input file (gene expression data you wish to apply STAMP on, in CSV format);
	   Name this file: ”gene_expression_input.csv”
            I. Make sure gene names are in Hugo symbol format, and are in a column named 'GENE_SYMBOLS’.
            II. Make sure samples are given as columns (The script will transpose this to generate predictions), 
	    	with the first row being the sample ids/names.
            III. These sample ids/names will be used in the output file to identify predictions of each sample.
        d. ‘models’ directory downloaded from our git, with the scripts it contains.

    4. To allow pytorch to work with sparse tensors, follow these instructions:
        a. Download the folder “Sparse_Tensors_Modified_Files/”.
        b. Go to the ‘torch/’ library in your conda environment. It should be at: anaconda3/envs/STAMP_env/lib/python3.7/site-packages/torch/.
        c. Modify the files torch/_utils.py, torch/tensor.py, torch.serialization.py: replace them with the files downloaded from our folder.
    	• These scripts’ modifications only add new functions to these torch scripts, allowing the saving and loading of models with sparse tensors 
		(as the GCN models are). 
	*** It will not hurt any other torch functionality. ***

    5. Run in command line:
	python Run_model_predictions_by_user_input.py


