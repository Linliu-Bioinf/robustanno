#!/bin/env python3

import os
import pandas as pd
import anndata
from scipy import sparse
import gzip
import numpy as np

"""
__date__ == "2022_05_11"
__author__ == "liulin4"
"""


def binExp(geneExpFile, binSize = 20):
    suffix = os.path.basename(geneExpFile).split('.')[-1]
    columnTypes = {"geneID": 'category', "x": np.uint32, "y": np.uint32, "values": np.uint32, "MIDCount": np.uint32, "UMICounts": np.uint32}
    if suffix == 'txt':
        df = pd.read_csv(geneExpFile, sep = '\t', dtype = columnTypes)
    elif suffix == 'pickle':
        df = pd.read_pickle(geneExpFile)
    elif suffix == 'tsv':
        df = pd.read_csv(geneExpFile, encoding = 'gbk', sep = '\t')
    elif suffix == "gz":
        header = parse_head(geneExpFile)
        df = pd.read_csv(geneExpFile, header = header, sep='\t', dtype=columnTypes)

    if 'values' in df.columns:
        df.rename(columns = {'values': 'MIDCounts'}, inplace = True)
    elif 'MIDCount' in df.columns:
        df.rename(columns = {'MIDCount': 'MIDCounts'}, inplace = True)
    elif 'UMICount' in df.columns:
        df.rename(columns = {'UMICount': 'MIDCounts'}, inplace = True)

    print("Start to do data binning!")
    df['x'] = (df['x']/int(binSize)).astype(np.uint32)
    df['y'] = (df['y']/int(binSize)).astype(np.uint32)
    df['binID'] = (df['x'].astype(str) + '-' + df['y'].astype(str)).astype('category')
    cell_list = df["binID"].astype('category')
    gene_list = df["geneID"].astype('category')
    data =  df["MIDCounts"].to_numpy()
    row = cell_list.cat.codes.to_numpy()
    col = gene_list.cat.codes.to_numpy()
    obs = pd.DataFrame(index = cell_list.cat.categories)
    var = pd.DataFrame(index = gene_list.cat.categories)
    X = sparse.csr_matrix((data, (row, col)), shape = (len(obs), len(var)))
    adata = anndata.AnnData(X, obs = obs, var = var)
    adata.obsm['spatial'] = pd.Series(adata.obs.index).str.split('-', expand=True).astype('int').values
    adata.raw = adata
    return adata

def parse_head(gem):
    if gem.endswith('.gz'):
        f = gzip.open(gem, 'rb')
    else:
        f = open(gem, 'rb')
    header = ''
    num_of_header_lines = 0
    eoh = 0
    for i, l in enumerate(f):
        l = l.decode("utf-8") 
        if l.startswith('#'):
            header += l
            num_of_header_lines += 1
            eoh = f.tell()
        else:
            break
    f.seek(eoh)
    return num_of_header_lines


