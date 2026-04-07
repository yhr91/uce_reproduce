# NOTE: this script makes running scGPT (zero-shot) a lot easier
# It is modified from https://github.com/microsoft/zero-shot-scfoundation
# NOTE: scGPT code for running highly variable genes can cause errors, we reccomend running highly variably genes separately, saving the anndata file (with the .X as dense)
# and then running this file with the saved path
# Make sure to update this script with the correct values for gene_col, batch_col, label_cols and layer_key based on your saved anndata!




import os
import logging
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

from sc_foundation_evals import cell_embeddings, scgpt_forward, data, model_output
from sc_foundation_evals.helpers.custom_logging import log

import argparse

# specify the path to anndata object
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='scGPT eval an anndata')
    parser.add_argument('--path', default="")
       
    
    args = parser.parse_args()
    dataset_path = args.path
    log.setLevel(logging.INFO)
    
    # path to the pre-trained model, 3 files are expected: 
    # model_weights (best_model.pt), model args (args.json), and model vocab (vocab.json)
    model_dir="../data/weights/scgpt/scGPT_human"
    # batch_size depends on available GPU memory; should be a multiple of 8
    batch_size=32
    # output_dir is the path to which the results should be saved
    output_dir="../output/scgpt/scgpt_human/"
    # path to where we will store the embeddings and other evaluation outputs
    model_out = os.path.join(output_dir, "model_outputs")
    # if you can use multithreading specify num_workers
    num_workers=0
    
    # path to the pre-trained model, 3 files are expected: 
    # model_weights (best_model.pt), model args (args.json), and model vocab (vocab.json)
    model_dir="../data/weights/scgpt/scGPT_human"
    # batch_size depends on available GPU memory; should be a multiple of 8
    batch_size=32
    # output_dir is the path to which the results should be saved
    output_dir="../output/scgpt/scgpt_human/"
    # path to where we will store the embeddings and other evaluation outputs
    model_out = os.path.join(output_dir, "model_outputs")
    # if you can use multithreading specify num_workers
    num_workers=0
    
    input_bins=51
    model_run="pretrained"
    seed=7
    n_hvg=3000
    # maximum sequence of the input is controlled by max_seq_len, here I'm using the pretrained default
    max_seq_len=3000
    
    
    
    
    # in which column in adata.obs are gene names stored? if they are in index, the index will be copied to a column with this name
    gene_col = "feature_name"
    # batch column found in adata.obs
    batch_col = "donor" #"batch"
    # where are labels stored in adata.obs? 
    label_cols = ["cell_ontology_class"]
    # where the raw counts are stored?
    layer_key = "counts"
    # are the values log_norm already?
    log_norm = False
    
    
    
    
    # create the model
    scgpt_model = scgpt_forward.scGPT_instance(saved_model_path = model_dir,
                                               model_run = model_run,
                                               batch_size = batch_size, 
                                               save_dir = output_dir,
                                               num_workers = num_workers, 
                                               explicit_save_dir = True)
    
    
    # create config
    scgpt_model.create_configs(seed = seed, 
                               max_seq_len = max_seq_len, 
                               n_bins = input_bins)
    
    
    scgpt_model.load_pretrained_model()
    
    
    input_data = data.InputData(adata_dataset_path = dataset_path)
    
    vocab_list = scgpt_model.vocab.get_stoi().keys()
    
    input_data.preprocess_data(gene_vocab = vocab_list,
                               model_type = "scGPT",
                               gene_col = gene_col,
                               data_is_raw = not log_norm,
                               counts_layer = layer_key, 
                               n_bins = input_bins,
                               fract_matching=0.3, # YANAY CHANGE
                               n_hvg = n_hvg)
    
    scgpt_model.tokenize_data(data = input_data,
                              input_layer_key = "X_binned",
                              include_zero_genes = False)
    
    scgpt_model.extract_embeddings(data = input_data)
    
    eval_ce = cell_embeddings.CellEmbeddingsEval(scgpt_model,
                                                 data = input_data,
                                                 output_dir = model_out,
                                                 label_key = label_cols,
                                                 batch_key = batch_col)
    
    eval_ce.data.adata.write(dataset_path.replace(".h5ad", "_scgpt.h5ad"))
#results = eval_ce.evaluate(n_cells = 1000)
#print(results)
