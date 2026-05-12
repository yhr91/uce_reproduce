### Universal Cell Embedding (Reproduce Results from Paper)

Download the `export_data` directory files [here](https://drive.google.com/file/d/1O9a0UkqPHT_LOXxf695c1J-SyHL47eGI/view?usp=sharing).

The `notebooks/figures` directory contains files specifically for figure 3 analyses. 


`evaluate_geneformer.py` and `evaluate_scgpt.py` contain code for running Geneformer and scGPT more easily in a zero shot setting, adapted from https://github.com/microsoft/zero-shot-scfoundation

For benchmarks, such as for tabula sapiens, in the `scib.ipynb` the correctly loaded file from `export_data` has the precomputed embeddings for all models benchmarked already.


## File Manifest

| File Name | Description |
| --------- | ----------- |
| figures/ | Directory contains code for figure 3 analysis. (see that directory for list of files) |
|benchmark_uce_saturn_samap.ipynb | Run the labal transfer benchmarking for cross species comparison of UCE vs SATURN vs SAMAP for the 4 new species datasets |
| check_cosine_similarity.ipynb | Check consistency of embeddings for UCE models for the same cell run with different seeds, and different cell types, using cosin similarity. Also includes the comparison for scGPT using different Differential Expression settings. |
| epo_case_study.ipynb | Analysis for figure 4 for the EPO Norn cell case study |
| epo_heatmaps.ipynb | Code to make the heatmaps for figure 4 EPO Norn cell case study analysis. |
| evaluate_geneformer.py | This script is used to create geneformer zero-shot embeddings for a anndata dataset | 
| evaluate_scgpt.py | This script is used to create scGPT zero-shot embeddings for a anndata dataset |
| fake_cells.ipynb | This notebook contains the analysis for the results of embedding cells with random gene expression |
| fly_dendrogram.ipynb | Notebook for creting the fly cell atlas dendrogram |
| fly_embedding_check.ipynb | Notebook for comparing model cell typing with logistic classifier for the fly cell atlas with UCE, versus random embeddings. |
| fly_enterocytes.ipynb | Notebook contains the analysis code for nearest neighbors of human gut cells and fly gut cells (enterocyte case study) |
| plot_loss_vs_distance.ipynb | Notebook for creating the evolutionary distance versus loss figure. |
| plot_training_curves.ipynb | Notebook for plotting the various training model ablation loss curve figures, and the tokenization ablation figure. |
| scgraph_benchmark.ipynb | Notebook contains the code to run the scGraph benchmark on tabula sapiens and Human Brain Cell Atlas |
| tabula_scArches.ipynb | Notebook contains the code for creating scArches embeddings of tabula sapiens. |
| tabula_scvi.ipynb | Notebook contains the code for creating scVI embeddings of tabula sapiens. |
