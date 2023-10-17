### General project utils

#### Load packages
import numpy as np
import pandas as pd
import scanpy as sc
import diffxpy.api as de
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt


import clustering_errors_utils
import AUC_scoring_utils



#### Define utils

## Function: k_means clustering
def compute_kmeans(adata, k):
    
    # Get the data matrix from the adata object
    data = adata.X
    
    # Perform k-means clustering
    kmeans = KMeans(n_clusters=k, random_state=0).fit(data)
    
    # Get the cluster labels and convert them to string data type
    labels_str = kmeans.labels_.astype(str)
    
    # Add the cluster labels to the adata object as a categorical column
    adata.obs['k_means'] = labels_str
    
    # Return the modified adata object
    return adata


# Function: ARI computation
#     """
#     Compute the Adjusted Rand Index (ARI) between two sets of cluster assignments contained
#     in the `.obs` attribute of AnnData objects.
#     
#     Parameters
#     ----------
#     adata1: AnnData object
#         The first AnnData object containing the first set of cluster assignments.
#     adata2: AnnData object
#         The second AnnData object containing the second set of cluster assignments.
#     clust_method: str
#         The key in the `.obs` attribute of the AnnData object containing the cluster assignments.
#         
#     Returns
#     -------
#     float
#         The Adjusted Rand Index between the two sets of cluster assignments.
#     """
def compute_ari(adata1, adata2, clust_method):

    # Extract cluster assignments from the two AnnData objects
    cluster_labels_1 = adata1.obs[clust_method]
    cluster_labels_2 = adata2.obs[clust_method]
    
    # Compute and return the Adjusted Rand Index
    return adjusted_rand_score(cluster_labels_1, cluster_labels_2)
    




def get_cluster_genes(adata, clust_method, test, clust_label=None):
    """
    This function retrieves the set of significantly differentially expressed gene names that are associated with each cluster in a given dataset
    using a clustering method and test specified by the user.

    Parameters
    ----------
    - adata: Anndata object - the dataset containing the gene expression data.
    - clust_method: str - the clustering method used to cluster the cells in the dataset.
    - test: str- the DE test method to be used. One of 't-test','t-test_overestim_var','wilcoxon' or 'diffxpy'.
    - clust_label : str - optional. If provided, DE testing is performed only for the specified cluster.

    Returns:
    - A named list of unique gene names that are associated with all clusters in the dataset.
    """

    # Initialize an empty dictionary to store the gene sets for each cluster
    gene_dict = {}

    # Get unique clusters
    clusters = pd.Series(adata.obs[clust_method]).unique()

    # Check if ' ' category is present in adata and remove it
    if 'unassigned                   ' in adata.obs[clust_method].values:
        adata = adata[adata.obs[clust_method] != 'unassigned                   ']

    # Loop through each cluster and get the top genes
    for cluster in clusters:

        # If clust_label is specified, only do DE testing: clust_label vs rest
        if clust_label:
            if cluster != clust_label:
                continue

        # Check if the current cluster contains only one sample
        if len(adata.obs[adata.obs[clust_method] == cluster]) <= 1:
            continue

        if test == 'diffxpy':
            # Create a column called groupings containing 1s for the current cluster and 0s for all others
            adata.obs['groupings'] = [1 if c == cluster else 0 for c in adata.obs[clust_method]]

            # Perform differential expression analysis using diffxpy for the current cluster
            print(f"Performing differential expression test for {cluster}")

            test_DE = de.test.two_sample(adata, grouping='groupings')
            result = test_DE.summary()

            # Sort the genes based on the absolute value of the log2 fold change
            result_df = pd.DataFrame(result.sort_values(by='log2fc', ascending=False))

            # Get the gene names from the result and add them to the set
            gene_names = set(result_df.loc[result_df['qval'] < 0.01, 'gene'])

            gene_dict[cluster] = gene_names

        else:
        	# Perform differential expression analysis using diffxpy for the current cluster
            print(f"Performing differential expression test for {cluster}")
            
            # Run the rank_genes_groups_df function for the current cluster
            sc.tl.rank_genes_groups(adata, clust_method, method=test, groups=[cluster])
            result = sc.get.rank_genes_groups_df(adata, group=cluster)

            # Filter genes based on adjusted p-values
            gene_names = set(result.loc[result['pvals_adj'] < 0.01, 'names'].values)

            gene_dict[cluster] = gene_names

    return gene_dict


# # Function: compute Jaccard similarity between genes sets
# def jaccard_similarity(set1, set2):
# # 	"""
# #     Compute Jaccard similarity between two sets of genes
# #     
# #     Parameters
# #     ----------
# #     - set1 (set): first set of genes
# #     - set2 (set): second set of genes
# #     
# #     Returns
# #     ----------
# #     - float: Jaccard similarity index between set1 and set2
# #     """
#     intersection = len(set1.intersection(set2))
#     union = len(set1.union(set2))
#     return intersection / union

def jaccard_similarity(dict1, dict2):
    # Initialize variables to store intersection and union lengths
    intersect_lengths = []
    union_lengths = []
    
    # Iterate over the named sets in dict1 and dict2 and compute Jaccard similarity for each matching set
    for cluster_name in set(dict1.keys()) & set(dict2.keys()):
        set1 = set(dict1[cluster_name])
        set2 = set(dict2[cluster_name])
        intersect_len = len(set1.intersection(set2))
        union_len = len(set1.union(set2))
        intersect_lengths.append(intersect_len)
        union_lengths.append(union_len)
    
    # Compute and return the average Jaccard similarity index across all matching sets
    return (np.mean(np.array(intersect_lengths) /np.array(union_lengths)))    
    
# Function: produce plot of ARI values vs. Jaccard similarity indices
def plot_ari_jaccard(ari_values, jaccard_indices, title, xlabel,ylabel):
#     """
#     This function plots a line graph of ARI values and corresponding Jaccard indices.
# 
#     Parameters
#     ----------
#     - ari_values: numpy array - an array of ARI values.
#     - jaccard_indices: numpy array - an array of corresponding Jaccard indices.
#     - title: str - the title of the plot.
# 
#     Returns
#     -------
#     None.
#     """
    plt.figure(figsize=(8, 6))
    plt.plot(ari_values, jaccard_indices, '-o')
    plt.xlabel('ARI')
    plt.ylabel('Jaccard Index')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()



# Function: main 
def main_test_DE_method(adata, percentages, clust_method, test_method, clust_label = None):
#     """
#     Main function to test the performance of a DE analysis method under different levels of clustering errors.
# 
#     Parameters
#     ----------
#     adata (AnnData object): Annotated data matrix.
#     percentages (list of floats): List of percentages of cells to assign to a different cluster.
#     clust_method (str): Clustering method to be used ('louvain' or 'leiden').
#     test_method (str): DE gene extraction method to be used ('t-test','t-test_overestim_var','wilcoxon','diffxpy').
#     clust_label (optional): if specified, the pipeline is only applied to the selected cluster.
#     
#     Plots
#     ----------
#     ARI values vs. Jaccard indices for the datasets with introduced clustering errors.
# 
#     Returns
#     ----------
#     Jaccard_indices (list of floats): List of Jaccard similarity indices between the set of DE genes extracted from the original
#     clustering and the sets of DE genes extracted from the datasets with introduced clustering errors. 
#     """
    
    # Introduce clustering errors
    clust_error_adata = []
    
    for p in percentages:
        print(f"Running clust_error with percentage={p}")
        # ce = clust_error(adata, p, clust_method, clust_label)
        ce = clust_error(adata, p, clust_method)
        clust_error_adata.append(ce)
    
    # Compute ARI and Extract DE_genes for each adata object in clust_error_adata
    ARIs = []
    DEGs = []
    
    for ce_adata in clust_error_adata:
        #print(f"computing ARI and extracting DEGs for ={ce_adata}")
        ARI_val = compute_ari(adata, ce_adata, clust_method)
        ARIs.append(ARI_val)
        
        DEGs_ce_adata_set = get_cluster_genes(ce_adata, clust_method, test_method, clust_label)
        DEGs.append(DEGs_ce_adata_set)
    
    # Compute Jaccard similarity index between original set of DEGs and all sets of DEGs in DEGs
    Jaccard_indices = []
    
    for DEGs_set in DEGs:
        similarity = jaccard_similarity(DEGs[0],DEGs_set)
        Jaccard_indices.append (similarity)
        
    # Produce final ARI vs. Jaccard plot
    if clust_label is not None:
        title = "ARI vs. Jaccard - " + test_method + ' ' + clust_label
        xlabel = "global ARI"
        ylabel = "cluster-specific Jaccard index"
    else:
        title = "ARI vs. Jaccard - " + test_method
        xlabel = "ARI"
        ylabel = "Jaccard index"

    plot_ari_jaccard(ARIs, Jaccard_indices, title, xlabel, ylabel)
  
    return ARIs,Jaccard_indices
    
### Pipeline wrappers: global vs single clust

def plot_ari_jaccard_methods_comparison(ari_values, jaccard_indices, comparison_names, title, xlabel, ylabel):
    """
    This function plots a line graph of average ARI values and corresponding Jaccard indices for each method.

    Parameters
    ----------
    - ari_values: list - a list of lists containing ARI values for each method.
    - jaccard_indices: list - a list of lists containing Jaccard indices for each method.
    - comparison_names: list - a list of instances to compare.
    - title: str - the title of the plot.
    - xlabel: str - the label for the x-axis.
    - ylabel: str - the label for the y-axis.

    Returns
    -------
    None.
    """
    plt.figure(figsize=(8, 6))
    
    # Create a color map for each method
    color_map = plt.cm.get_cmap('tab10')
    
    # Plot data for each method
    for i, method in enumerate(comparison_names):
        ari_vals = ari_values[i]
        jaccard_vals = jaccard_indices[i]
        
        # Plot ARI vs. Jaccard for the current method
        plt.plot(ari_vals, jaccard_vals, '-o', label=method, color=color_map(i / len(comparison_names)))
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.show()


def plot_ari_jaccard_methods_comparison_wrapper(adata, percentages, clust_method, test_methods, clust_label=None):
    """
    This function calls the `main_test_DE_method` for each method in the input list and collects the ARIs and Jaccard indices.
    Then it sends them as input to the `plot_ari_jaccard_methods_comparison` function to obtain the final plot.

    Parameters
    ----------
    - adata: AnnData object - Annotated data matrix.
    - percentages: list of floats - List of percentages of cells to assign to a different cluster.
    - clust_method: str - Clustering method to be used ('louvain' or 'leiden').
    - test_methods: list of str - DE gene extraction methods to be used.
    - clust_label: optional - If specified, the pipeline is only applied to the selected cluster.

    Returns
    -------
    None.
    """
    ARIs = []
    Jaccard_indices = []
    comparison_names = []
    
    for method in test_methods:
        print(f"Running main_test_DE_method for method={method}")
        ari_values, jaccard_values = main_test_DE_method(adata, percentages, clust_method, method, clust_label)
        ARIs.append(ari_values)
        Jaccard_indices.append(jaccard_values)
        comparison_names.append(method)

    title = "ARI vs. Jaccard - Average methods comparison"
    xlabel = "ARI"
    ylabel = "Jaccard index"

    plot_ari_jaccard_methods_comparison(ARIs, Jaccard_indices, comparison_names, title, xlabel, ylabel)


def plot_ari_jaccard_auto_clust_labels_comparison_wrapper(adata, percentages, clust_method, test_method):
    """
    This function retrieves the cluster labels from `adata.obs[clust_method]`, calls the `main_test_DE_method`
    for each cluster label, and collects the ARIs and Jaccard indices. Then it sends them as input to the
    `plot_ari_jaccard_methods_comparison` function to obtain the final plot.

    Parameters
    ----------
    - adata: AnnData object - Annotated data matrix.
    - percentages: list of floats - List of percentages of cells to assign to a different cluster.
    - clust_method: str - Clustering method to be used ('louvain' or 'leiden').
    - test_method: str - DE gene extraction method to be used.

    Returns
    -------
    None.
    """
    clust_labels = adata.obs[clust_method].unique().tolist()
    
    ARIs = []
    Jaccard_indices = []
    comparison_names = []
    
    for clust_label in clust_labels:
        print(f"Running main_test_DE_method for clust_label={clust_label}")
        ari_values, jaccard_values = main_test_DE_method(adata, percentages, clust_method, test_method, clust_label)
        ARIs.append(ari_values)
        Jaccard_indices.append(jaccard_values)
        comparison_names.append(test_method + " " + str(clust_label))
    
    title = "ARI vs. Jaccard - Cluster-specific comparison"
    xlabel = "Global ARI"
    ylabel = "Cluster-specific Jaccard index"

    plot_ari_jaccard_methods_comparison(ARIs, Jaccard_indices, comparison_names, title, xlabel, ylabel)