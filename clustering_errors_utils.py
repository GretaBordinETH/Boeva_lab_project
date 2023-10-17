### Clustering errors-introduction utils

#### Load packages
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial.distance import cdist


#### Define utils

## Function: introduction of clustering errors
def compute_centroids(adata, clusters, clust_method):
    """
    Compute the centroids for each cluster based on the PCA projection of the data.
    
    Parameters:
    ----------
        adata (AnnData): AnnData object containing the gene expression data.
        clusters (list): List of cluster labels.
        clust_method (str): Clustering method used to assign cell labels to the AnnData object.
    
    Returns:
    ----------
        centroids (dict): Dictionary containing the centroid coordinates for each cluster label.
    """
    
    # initialize empty centroids list
    centroids = {}
    pca_df = pd.DataFrame(adata.obsm['X_pca'][:, :2], index=adata.obs.index)
    
    # calculate centroid for each cluster
    for cluster in clusters:
        # subset the data to get only the points in the cluster
        cluster_points = adata.obs.index[adata.obs[clust_method] == cluster]
        df_cluster = pca_df.loc[cluster_points]

        # compute the centroid
        centroid = df_cluster.mean(axis=0)
        centroids[cluster] = centroid

    return centroids
    
## Function: compute the distance of each cell from other clusters centroids
def compute_min_centroids_dists(adata, clust_method):
	"""
    Compute the minimum distance between each cell in a dataset and its closest cluster centroid, return the
    identity of the closest cluster for each cell as well as its distance from it.

    Parameters:
    ----------
        - adata (AnnData object): Annotated data matrix containing the gene expression data and cell metadata.
        - clust_method (str): Name of the column in adata.obs that contains the cluster labels.

    Returns:
    ----------
        - A dictionary with two keys:
            - "min_dists": A 1D numpy array of shape (n_cells,) containing the minimum distance between each cell and its closest centroid.
            - "closest_cluster": A 1D numpy array of shape (n_cells,) containing the identity of the closest cluster for each cell.
    """

    # get list of clusters
	clusters = sorted(adata.obs[clust_method].unique())

    # get df of dim_red coordinates
	pca_df = pd.DataFrame(adata.obsm['X_pca'][:, :2], index=adata.obs.index)

    # compute cluster centroids
	centroids = compute_centroids(adata, clusters, clust_method)

    # initialize output vectors
	min_dists = np.repeat(np.Inf, pca_df.shape[0])
	closest_cluster = np.repeat('unassigned', pca_df.shape[0])

    # loop over cells and compute distances to other centroids
	for i, cell in enumerate(pca_df.values):
        # get cluster label for current cell
		cell_cluster = adata.obs[clust_method][i]

        # loop over all other centroids and compute distances
		for j in [c for c in clusters if c != cell_cluster]:
			dist = cdist(cell.reshape(1, -1), centroids[j].values.reshape(1, -1))[0][0]
			if dist < min_dists[i]:
				min_dists[i] = dist
				closest_cluster[i] = j

	return {"min_dists": min_dists, "closest_cluster": closest_cluster}
	

## Function: compute average connectivity to cells in other clusters
def compute_max_connectivity(adata, clust_method):
    """
    Compute the max average connectivity of each cell to the cells in the other clusters, and the corresponding 
    cluster label. The connectivity information is stored in adata.obsp['connectivities'].

    Parameters:
    ----------
        - adata (AnnData object): Annotated data matrix containing the gene expression data and cell metadata.
        - clust_method (str): Name of the column in adata.obs that contains the cluster labels.

    Returns:
    ----------
        - A dictionary with two keys:
            - "max_connectivity": A 1D numpy array of shape (n_cells,) containing the max average connectivity of each cell to the cells in the other clusters.
            - "most_connected_cluster": A 1D numpy array of shape (n_cells,) containing the cluster label of the cluster that the cell is most connected to.
    """
    
    # get list of clusters
    clusters = sorted(adata.obs[clust_method].unique())

    # initialize output vectors
    max_connectivity = np.repeat(-1.0, adata.shape[0])
    most_connected_cluster = np.repeat('unassigned', adata.shape[0])

    # loop over cells and compute max average connectivity to cells in other clusters
    for i, cell in enumerate(adata.obs.index):
        # get cluster label for current cell
        cell_cluster = adata.obs[clust_method][i]
        
        # initialize dictionary to store connectivity to cells in each cluster
        cluster_connectivity = {c: [] for c in clusters}
        
        # loop over neighbors of current cell and compute connectivity
        neighbors = adata.obsp['connectivities'][i].nonzero()[1]
        for neighbor in neighbors:
            neighbor_cluster = adata.obs[clust_method][neighbor]
            if neighbor_cluster != cell_cluster:
                cluster_connectivity[neighbor_cluster].append(adata.obsp['connectivities'][i, neighbor])
        
        # compute max average connectivity and corresponding cluster label
        for c in clusters:
            if len(cluster_connectivity[c]) > 0:
                avg_connectivity = np.mean(cluster_connectivity[c])
                if avg_connectivity > max_connectivity[i]:
                    max_connectivity[i] = avg_connectivity
                    most_connected_cluster[i] = c
    
    # return output as a dictionary
    return {"max_connectivity": max_connectivity, "most_connected_cluster": most_connected_cluster}



## Function: introduce clustering errors
def clust_error(adata, percent, clust_method):#, clust_label=None):
# 	"""
# 	Switches a given percentage of cell cluster assignments in adata and returns the resulting AnnData object.
# 
# 	Parameters
# 	----------
# 	adata : AnnData object
# 		The annotated data matrix of shape (n_obs, n_vars).
# 	percent : float
# 		The percentage of cells whose cluster assignment should be switched. Must be a value between 0 and 100.
# 	clust_method : str
# 		The name of the column in adata.obs that contains the cell cluster assignments.
# 	# clust_label : str - optional
# 	If provided, clustering errors are only introduced in the specified cluster.
# 
# 	Returns
# 	-------
# 	error_adata : AnnData object
# 		The modified AnnData object with switched cell cluster assignments.
# 	"""

    error_adata = adata.copy()
    
    n_cells = adata.shape[0] # if not clust_label else sum(adata.obs[clust_method] == clust_label)
    n_cells_to_switch = int((n_cells * percent) / 100)
    
    
    ### KNN connectivity-based cluster perturbation
    
    if clust_method in ['leiden', 'louvain']:
        max_connectivities, most_connected_clusters = compute_max_connectivity(adata, clust_method).values()
        error_adata.obs['most_connected_cluster'] = most_connected_clusters
        error_adata.obs['max_connectivity'] = max_connectivities
 
        # sort cells by max connectivity
        sorted_cells = np.argsort(max_connectivities)[::-1]
        
        # select cells for which cluster assignments have to be switched
        cells_to_switch = sorted_cells[:n_cells_to_switch]
        
#         if clust_label:
#             clust_indices = np.where(error_adata.obs[clust_method] == clust_label)[0]
#             cells_to_switch = clust_indices[:n_cells_to_switch]
#         else:
#             cells_to_switch = sorted_cells[:n_cells_to_switch]
    

   ### closest centroid in PCA space cluster perturbation
    else:
        min_dists, closest_clusters = compute_min_centroids_dists(adata, clust_method).values()
        error_adata.obs['dist'] = min_dists
        error_adata.obs['closest_centroid'] = closest_clusters
        
		# sort cells by min distance to closest centroid
        sorted_cells = np.argsort(min_dists)
        
        # select cells for which cluster assignments have to be switched
        cells_to_switch = sorted_cells[:n_cells_to_switch]
        
#         if clust_label:
#             clust_indices = np.where(error_adata.obs[clust_method] == clust_label)[0]
#             cells_to_switch = clust_indices[:n_cells_to_switch]
#         else:
#             cells_to_switch = sorted_cells[:n_cells_to_switch]
    
    for i in cells_to_switch:
        if clust_method in ['leiden', 'louvain']:
            new_cluster = error_adata.obs['most_connected_cluster'][i]
        else:
            new_cluster = error_adata.obs['closest_centroid'][i]
        
        if new_cluster not in error_adata.obs[clust_method].cat.categories:
            error_adata.obs[clust_method] = error_adata.obs[clust_method].astype(str)
            error_adata.obs[clust_method][i] = str(new_cluster)
            error_adata.obs[clust_method] = pd.Categorical(error_adata.obs[clust_method])
        else:
            error_adata.obs[clust_method][i] = new_cluster
            
    return error_adata