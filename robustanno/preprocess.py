#!/bin/env python3

import scanpy as sc
import numpy as np
import scipy

"""
__date__ ==  2022-05-19
__author__ == "liulin4"
"""

def preprocess(adata, pre = "CLR", dataType = "SP"):
    adata.var_names_make_unique()
    if dataType == "SP":
        adata.obs["total_counts"] = np.ravel(adata.X.sum(axis = 1)) 
    else:
        adata.var["mt"] = adata.var_names.str.startswith(("mt-", "Mt-", "MT-"))
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
        adata = adata[adata.obs["pct_counts_mt"] < 10]
    
    max_counts = np.percentile(np.array(adata.X.sum(axis = 1)).ravel(), 98)
    min_counts = np.percentile(np.array(adata.X.sum(axis = 1)).ravel(), 2)
    sc.pp.filter_cells(adata, max_counts = max_counts, inplace = True)
    sc.pp.filter_cells(adata, min_counts = min_counts, inplace = True)
    sc.pp.filter_genes(adata, min_cells = 10)
    if pre == "CLR":
        adata = clr_normalization(adata)
        print("Centered log-ratio normalization is done")
    elif pre == "ST":
        sc.pp.normalize_per_cell(adata, counts_per_cell_after = adata.obs["total_counts"].median())
        sc.pp.log1p(adata)
        sc.pp.scale(adata, max_value = 10)
    else:
        print("No normalization is done!")

    return adata


def clr_normalization(adata):
    #Centered Log-Ratio normalization (CLR): clr(i) = ln(x_i/g(x)); g(x) = n^(x_2 ...x_n)  
    if scipy.sparse.issparse(adata.X):
        GM = np.exp(np.sum(np.log(adata.X.toarray() + 1)/adata.X.shape[1], axis = 1))
        adata.X = np.log(adata.X.toarray()/GM[:, None] + 1)
    else:
        GM = np.exp(np.sum(np.log(adata.X + 1)/adata.X.shape[1], axis = 1))
        adata.X = np.log(adata.X/GM[:, None] + 1)
    return adata
