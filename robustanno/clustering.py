#!/bin/env python3

import os
import numpy as np
import scanpy as sc
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

"""
__date__ == 2022-05-13
__author__ == liulin4
"""

def clustering(adata, method = "leiden", n_clusters = None, n_neighbors = 15, n_pcs = 30, r = 0.2, num_filter = 20):
    sc.tl.pca(adata, svd_solver = "arpack")
    sc.pp.neighbors(adata, n_neighbors = n_neighbors, n_pcs = n_pcs)

    if method == "kmeans":
        assert n_clusters > 1
        from sklearn.cluster import KMeans
        kmeans = KMeans(n_clusters = n_clusters, random_state = 0).fit(adata.obsm['X_pca'][:,:n_pcs])
        adata.obs[method] = kmeans.labels_.astype(str)
        print("{} clusters are saved in kmeans object!".format(n_clusters))
    else:
        print("Clustering by leiden!")
        sc.tl.leiden(adata, key_added = "leiden", resolution = r)
    sc.tl.umap(adata)

    adata = adata[adata.obs[method].isin(pd.DataFrame(adata.obs[method].value_counts() >= num_filter).index.tolist())]
    print("Clusters less than {} cells are removed!".format(num_filter))
    return adata

