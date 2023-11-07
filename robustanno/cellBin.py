import numpy as np
import pandas as pd
import anndata
from scipy import sparse
import scanpy as sc
import os, sys
from collections import defaultdict

"""
__date__ == "2021-12-30"
__author__ == "liulin4@genomics.cn"
"""

def cellBin(geneExpFile):
    columnTypes = {"geneID": 'category', "x": np.uint32, "y": np.uint32, "UMICount": np.uint32, 'label': 'category', 'tag': 'category'}
    df = pd.read_csv(geneExpFile, sep='\t', dtype = columnTypes)

    cell_list = df['label'].astype('category')
    gene_list = df['geneID'].astype('category')

    data = df['UMICount'].to_numpy()
    row = cell_list.cat.codes.to_numpy()
    col = gene_list.cat.codes.to_numpy()
    obs = pd.DataFrame(index=cell_list.cat.categories)
    var = pd.DataFrame(index=gene_list.cat.categories)
    X = sparse.csr_matrix((data, (row, col)), shape=(len(obs), len(var)))
    adata = anndata.AnnData(X, obs=obs, var=var)
    cell_coo = df.groupby('label').mean()[['x','y']]
    adata.obsm['spatial'] = cell_coo.to_numpy()
    adata.raw = adata
    return adata


