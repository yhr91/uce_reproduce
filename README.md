### Universal Cell Embedding (Reproduce Results from Paper)

Download the `export_data` directory files [here](https://drive.google.com/file/d/1O9a0UkqPHT_LOXxf695c1J-SyHL47eGI/view?usp=sharing).

The `notebooks/figures` directory contains files specifically for figure 3 analyses. 


`evaluate_geneformer.py` and `evaluate_scgpt.py` contain code for running Geneformer and scGPT more easily in a zero shot setting, adapted from https://github.com/microsoft/zero-shot-scfoundation

For benchmarks, such as for tabula sapiens, in the `scib.ipynb` the correctly loaded file from `export_data` has the precomputed embeddings for all models benchmarked already.