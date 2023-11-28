import os
import glob
import pickle as pkl
import numpy as np
import scanpy as sc

def read_uce_embed(name, model=None):
    # Define the directory path
    if model == '4layer':
        dir_path = '/dfs/project/uce/model_output/uce_all_proc_v2/4layer/'
    elif model == '33layer':
        dir_path = '/lfs/local/0/yanay/uce_all_proc_33_8ep/'
    else:
        raise DataError("Wrong argument")
    
    # Use glob to find matching files
    matching_files = glob.glob(dir_path + '*' + name + '*.npz')
    
    if not matching_files:
        raise FileNotFoundError(f"No matching .npz file found for '{name}' in '{dir_path}'")
    
    fname = matching_files[0]
    ncells = int(fname.split('/')[-1].split('_')[1])

    # Load UCE Data
    if not os.path.exists(fname):
        raise FileNotFoundError(f"File not found: {fname}")
    
    npz_path = fname
    arr = np.memmap(npz_path, shape=(ncells, 1280), dtype="float32", mode='r')
    arr = np.array(arr)

    # Load Cell Types
    ct_path = f"/lfs/local/0/yanay/all_cell_types/{name}.pkl"
    
    if not os.path.exists(ct_path):
        raise FileNotFoundError(f"File not found: {ct_path}")
    
    with open(ct_path, "rb") as f:
        cell_types = pkl.load(f)
        
    adata = sc.AnnData(arr)
    adata.obs['cell_type'] = cell_types
    return adata
