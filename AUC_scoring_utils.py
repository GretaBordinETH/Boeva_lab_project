### AUC scores from scanpy gene-scoring utils

#### Load packages
import numpy as np
import pandas as pd
import scanpy as sc
import sklearn.metrics

#### Define utils

### def get_cluster_gene_score function
def get_cluster_gene_score(adata, clust_label, gene_list):
 
    # Compute gene scores
    sc.tl.score_genes(adata, gene_list=gene_list)
    
    # Assign gene scores to a column in adata.obs
    adata.obs['score_' + clust_label] = adata.obs['score']
    
    return adata

### def get cluster_auc_score function
def get_cluster_auc_score(adata, clust_label):
    
    # Create column names for cluster membership and score based on the cluster label
    membership_col = clust_label + '_membership'
    score_col = 'score_' + clust_label
    
    # Compute AUC score
    auc_score = sklearn.metrics.roc_auc_score(adata.obs[membership_col], adata.obs[score_col])
    
    return auc_score

### def get avg AUC score per dataset
def get_AUC(adata, clust_method, de_genes_dict, clust_label = None):
    print('Retrieving avg AUC')
    
    # Iterate over unique cluster labels
    unique_clusters = np.unique(adata.obs[clust_method]) if not clust_label else [clust_label]
    
    
    # Initialise empty list of AUC scores
    auc_scores = []
    
    for clust_label in unique_clusters:
        
        # Get differentially expressed genes for the cluster from dictionary of DEGs per cluster
        gene_list = de_genes_dict.get(clust_label, [])
        
        # Checking if the current cluster has no DE genes -> not computing gene scores
        if len(gene_list)== 0:
            continue

        
        # Compute gene scores for the cluster and add it to the adata object
        adata = get_cluster_gene_score(adata, clust_label, gene_list)
        
        # Create membership column for the cluster and add it to the adata object
        membership_col = clust_label + '_membership'
        adata.obs[membership_col] = np.where(adata.obs[clust_method] == clust_label, 1, 0)
        
        # Check if the new column only contains 0s
        if np.all(adata.obs[membership_col] == 0):
            print("Warning: The", membership_col, "column contains only 0s.")
        
        # Compute AUC score for the cluster
        auc_score = get_cluster_auc_score(adata, clust_label)
        print(clust_label)
        print(auc_score)
        auc_scores.append(auc_score)
    
    # Calculate average AUC score over all clusters
    avg_auc_score = np.mean(auc_scores)
    
    return avg_auc_score