#!/bin/env python3

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from collections import Counter
import warnings
warnings.filterwarnings("ignore")

"""
__date__ == 2022-05-20
__author__ == liulin4
markers has to be provided!
"""

def anno(adata, marker_info, key = "leiden", n_genes = 3, rank_method = "ROC", key_add = "celltype_pred"):
    if type(marker_info) is dict:
        df_marker = convert_dict2df(marker_info)
    elif isinstance(marker_info, pd.DataFrame):
        df_marker = marker_info
    else:
        TypeError("marker_info has to be provided as a dataframe or dict format!")

    sc.tl.rank_genes_groups(adata, key, method = "wilcoxon")
#    plt.rcParams["figure.figsize"] = [6,6]
#    sc.pl.rank_genes_groups(adata, n_genes = n_genes, sharey=False,save = "leiden_deg.png")
    rank_mtx = np.zeros((adata.obs[key].nunique(), len(df_marker.columns)))
    for idx, item in enumerate(adata.obs[key].cat.categories):
        dep = adata.uns["rank_genes_groups"]['names'][item][:n_genes].tolist()
        for idx_, g in enumerate(df_marker.columns):
            if g in dep:
                rank_mtx[idx, idx_] = dep.index(g) + 1
                
    w = weighting(rank_mtx, N = n_genes, method = rank_method) #row: marker; col: cluster
    scores = np.matmul(df_marker.to_numpy(), w)
    pred = df_marker.index[np.argmax(scores, axis = 0)].tolist() 
    marker_idx = np.argmax(w, axis = 0)
    m = df_marker.columns[marker_idx].tolist()
    type_marker = list(map(": ".join, zip(pred, m)))
    print("Type and related decisive marker: {}".format(type_marker))    

    deps = [adata.uns["rank_genes_groups"]['names'][i][0] for i in adata.obs[key].cat.categories]
    pred_type = list(map(": ".join, zip(pred, deps))) 

    # For clusters do not have markers in deps
    unknown_cluster = np.where(~scores.any(axis=0))[0].tolist()
    if len(unknown_cluster) > 0:
        print("No markers are differently expressed in clusters: {}".format(unknown_cluster))
        for i in unknown_cluster:
            pred_type[i] = "Unknown: " + "_".join([i for i in adata.uns["rank_genes_groups"]['names'][str(i)][0:3]])
    # Add second marker for clusters annotated as the duplicated type
    pred_type = check_rename_dup(adata, pred_type)

    adata.obs[key_add] = adata.obs[key]
    for idx, item in enumerate(adata.obs[key].cat.categories):
        adata.obs[key_add].replace(item, pred_type[idx], inplace = True)
    return adata


def check_rename_dup(adata, pred_type):
    kth_dep = 1
    while len(pred_type) > len(set(pred_type)):
        dup_dict = dict(Counter(pred_type))
        dup_type = [t for t, v in dup_dict.items() if v > 1]
        for t in dup_type:
            dup_idx = [idx for (idx, val) in enumerate(pred_type) if val == t]
            for k in dup_idx:
                pred_type[k] = t + "_" + adata.uns["rank_genes_groups"]['names'][str(k)][kth_dep]
            kth_dep += 1
    else:
        print("Annotated types are unique!")
    return pred_type

    
def weighting(rank_mtx, N = 3, method = "ROC"):
    #rank_mtx: float
    #ROC: (1/N) * sum_{k=n}^N(1/r_k)
    #RR: (1/r_n) * (1/sum_{k=1}^N (1/r_k))
    weights = np.zeros_like(rank_mtx)
    n, m = rank_mtx.shape
    rank_mtx = rank_mtx.astype(int)
    if method == "RR":
        deno = 0
        for k in range(1, N+ 1):
            deno += 1/k     
        for i in range(n):
            for j in range(m):
                if rank_mtx[i, j] != 0:
                    weights[i, j] = 1/rank_mtx[i, j]/deno
    elif method == "ROC":
        deno = N
        for i in range(n):
            for j in range(m):
                if rank_mtx[i, j] != 0:
                    numerator = 0
                    for k in range(rank_mtx[i,j], N+1):
                        numerator += 1/k
                    weights[i, j] =  numerator/N
    return weights.T 


def convert_dict2df(marker_dict):
    markers = []
    for i in list(marker_dict.values()):
        for j in i:
            markers.append(j)   
    df_marker = pd.DataFrame(np.zeros((len(marker_dict.keys()), len(markers))), index = marker_dict.keys(), columns = markers).astype(int)
    for label, gene in marker_dict.items():
        df_marker.loc[label, gene] = 1
    return df_marker




def celltypeAnno_classic(adata, marker_dict, key = "leiden", n_genes = 3, key_add = "celltype_pred"):
    sc.tl.rank_genes_groups(adata, key, method = "wilcoxon")
    deg_dict = {}
    for i in adata.obs[key].cat.categories:
        deg_dict[i] = adata.uns["rank_genes_groups"]['names'][i][:n_genes].tolist()

    n_cluster = adata.obs[key].nunique()    
    n_type = len(marker_dict.keys())
    df_score = pd.DataFrame(np.zeros((n_cluster, n_type)), index = adata.obs[key].cat.categories, columns = marker_dict.keys())
    for label, marker in marker_dict.items():
        g_type = len(marker)
        s = []
        for j in adata.obs[key].cat.categories:
            l_type = len(set(marker_dict[label]).intersection(set(deg_dict[j])))
            s.append(l_type)
        s = np.array(s)
        tmpdf = pd.DataFrame(np.sqrt(s/(s+1) * g_type/(g_type+1)), index = adata.obs[key].cat.categories, columns = [label]) 
        df_score.update(tmpdf)
    df_anno = pd.DataFrame(df_score.columns[np.argmax(df_score.values, axis = 1)].values, index = adata.obs[key].cat.categories, columns = ["cellType"])
    df_anno["scoring"] = np.max(df_score.values, axis = 1)

    adata.obs[key_add] = adata.obs[key]
    for i in adata.obs[key].cat.categories:
        adata.obs[key_add].replace(i, df_anno.loc[i, "cellType"], inplace = True)
    return adata, df_anno



