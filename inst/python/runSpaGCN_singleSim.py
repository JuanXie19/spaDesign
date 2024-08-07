import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg
import louvain
import sklearn
import os
import argparse


def main(input_count, input_spatial, output_file):
    # read the input files
    adata = sc.read_csv(input_count, delimiter = '\t')
    spatial = pd.read_csv(input_spatial, delimiter = ',')
    
    # check the dimension of the input count matrix
    print("Dimensions of input count DataFrame:", adata.shape)
    
    # checking the column names for debugging
    print("Columns in spatial DataFrame:", spatial.columns)
    
    # processing the input
    adata.obs['cluster'] = list(spatial['regions'])
    adata.obs['x'] = list(spatial['x'])
    adata.obs['y'] = list(spatial['y'])
    adata.var_names=[i.upper() for i in list(adata.var_names)]
    adata.var["genename"]=adata.var.index.astype("str")
    
    #Set coordinates
    x_array=adata.obs["x"].tolist()
    y_array=adata.obs["y"].tolist()
    
    #Calculate adjacent matrix
    s=1
    b=49
    adj=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
    adata.var_names_make_unique()
    spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
    spg.prefilter_specialgenes(adata)
    
    #Normalize and take log for UMI
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    p=0.5 
    
    #Find the l value given p
    l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)
    
    #For this toy data, we set the number of clusters=6 since this tissue has 6 layers
    n_clusters=spatial['regions'].nunique()
    
    #Set seed
    r_seed=t_seed=n_seed=100
    
    #Search for suitable resolution
    res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)
    clf=spg.SpaGCN()
    clf.set_l(l)
    
    #Set seed
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    
    #Run
    clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
    y_pred, prob=clf.predict()
    adata.obs["pred"]= y_pred
    adata.obs["pred"]=adata.obs["pred"].astype('category')
    
    # calculate ARI
    #ARI = sklearn.metrics.adjusted_rand_score(adata.obs['cluster'],adata.obs['pred'])
    #ARI
    df = adata.obs[['x','y','cluster','pred']]
    df.to_csv(output_file)
    print(f"Processed data saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process spatial transcriptomics data with SpaGCN")
    parser.add_argument("input_count", type=str, help="Path to the input count file")
    parser.add_argument("input_spatial", type=str, help="Path to the input spatial location file")
    parser.add_argument("output_file", type=str, help="Path to save the output CSV file")
    
    args = parser.parse_args()
    main(args.input_count, args.input_spatial, args.output_file)
