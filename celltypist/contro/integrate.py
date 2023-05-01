import numpy as np
import pandas as pd
from anndata import AnnData
from typing import Optional, Union
import types
import sys
from sklearn.neighbors import DistanceMetric
from sklearn.neighbors import KDTree
from sklearn.metrics import pairwise_distances
from .. import logger
from .symbols import UNASSIGN
from umap.umap_ import fuzzy_simplicial_set
import pynndescent
from annoy import AnnoyIndex
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
try:
    import faiss
except ImportError:
    pass

def _locate_meta_neighbors(pca, batch_list, celltype_list, n_meta_neighbors):
    """
    Find the cell type neighbors for a given cell type (union across batches).
    """
    celltypes = np.unique(celltype_list)
    if n_meta_neighbors == 1 or len(celltypes) == 1:
        cns = {celltype: np.array([celltype], dtype = object) for celltype in celltypes}
    else:
        c2b = {celltype: np.unique(batch_list[celltype_list == celltype]) for celltype in celltypes}
        b2d = {}
        for batch in np.unique(batch_list):
            batch_flag = batch_list == batch
            sub_celltypes = np.unique(celltype_list[batch_flag])
            dists = pairwise_distances(np.array([pca[batch_flag & (celltype_list == sub_celltype)].mean(axis=0) for sub_celltype in sub_celltypes]))
            b2d[batch] = pd.DataFrame(dists, index = sub_celltypes, columns = sub_celltypes)
        cns = {celltype: np.unique(np.concatenate([b2d[b].loc[celltype].sort_values().index[:n_meta_neighbors] for b in c2b[celltype]])) for celltype in celltypes}
    if UNASSIGN in celltypes:
        cns_keys = list(cns.keys())
        for cns_key in cns_keys:
            if cns_key == UNASSIGN:
                cns[cns_key] = celltypes
            else:
                cns[cns_key] = np.unique(np.append(cns[cns_key], UNASSIGN))
    return cns

def _create_tree(pca, computation, metric, annoy_n_trees, pynndescent_n_neighbors, pynndescent_random_state):
    """
    Copied from BBKNN. Create a faiss/cKDTree/KDTree/annoy/pynndescent index for nearest neighbor lookup.
    """
    if computation == 'annoy':
        ckd = AnnoyIndex(pca.shape[1], metric = metric)
        for i in np.arange(pca.shape[0]):
            ckd.add_item(i, pca[i, :])
        ckd.build(annoy_n_trees)
    elif computation == 'pynndescent':
        ckd = pynndescent.NNDescent(pca, metric = metric, n_jobs = -1, n_neighbors = pynndescent_n_neighbors, random_state = pynndescent_random_state)
        ckd.prepare()
    elif computation == 'faiss':
        ckd = faiss.IndexFlatL2(pca.shape[1])
        ckd.add(pca)
    elif computation == 'cKDTree':
        ckd = cKDTree(pca)
    elif computation == 'KDTree':
        ckd = KDTree(pca, metric = metric)
    return ckd

def _query_tree(pca, computation, ckd, n_neighbors):
    """
    Copied from BBKNN. Query the faiss/cKDTree/KDTree/annoy/pynndescent index.
    """
    if computation == 'annoy':
        ckdo_ind = []
        ckdo_dist = []
        for i in np.arange(pca.shape[0]):
            holder = ckd.get_nns_by_vector(pca[i, :], n_neighbors, include_distances = True)
            ckdo_ind.append(holder[0])
            ckdo_dist.append(holder[1])
        ckdout = (np.asarray(ckdo_dist), np.asarray(ckdo_ind))
    elif computation == 'pynndescent':
        ckdout = ckd.query(pca, k = n_neighbors)
        ckdout = (ckdout[1], ckdout[0])
    elif computation == 'faiss':
        D, I = ckd.search(pca, n_neighbors)
        D[D < 0] = 0
        ckdout = (np.sqrt(D), I)
    elif computation == 'cKDTree':
        ckdout = ckd.query(x = pca, k = n_neighbors, workers = -1)
    elif computation == 'KDTree':
        ckdout = ckd.query(pca, k = n_neighbors)
    return ckdout

def _get_graph(pca, batch_list, celltype_list, computation, n_neighbors, n_meta_neighbors, metric, annoy_n_trees, pynndescent_n_neighbors, pynndescent_random_state, random_state):
    """
    Identify the cell-type-controlled KNN structure to be used in graph construction.
    """
    celltype_groups = _locate_meta_neighbors(pca, batch_list, celltype_list, n_meta_neighbors)
    if computation == 'faiss':
        pca = pca.astype('float32')
    knn_dists = np.zeros((pca.shape[0], n_neighbors))
    knn_indices = np.copy(knn_dists).astype(int)
    #main
    celltypes = np.unique(celltype_list)
    for celltype in celltypes:
        flag_celltype = celltype_list == celltype
        ind_celltype = np.arange(len(batch_list))[flag_celltype]
        flag_group = np.isin(celltype_list, celltype_groups[celltype])
        batches = np.unique(batch_list[flag_group])
        q, mod = divmod(n_neighbors, len(batches))
        np.random.seed(random_state)
        n_neighbors_across = np.random.permutation([q]*(len(batches)-mod) + [q+1]*mod)
        cum_sum = np.cumsum(n_neighbors_across)
        for rank, batch, n_neighbors_each in zip(range(len(batches)), batches, n_neighbors_across):
            flag_batch = batch_list == batch
            flag = flag_group & flag_batch
            if n_neighbors_each > flag.sum():
                flag = flag_batch
            ind = np.arange(len(batch_list))[flag]
            ckd = _create_tree(pca[flag], computation, metric, annoy_n_trees, pynndescent_n_neighbors, pynndescent_random_state)
            ckdout = _query_tree(pca[flag_celltype], computation, ckd, n_neighbors_each)
            col_range = np.arange(0 if rank == 0 else cum_sum[rank-1], cum_sum[rank])
            knn_indices[ind_celltype[:, np.newaxis], col_range] = ind[ckdout[1]]
            knn_dists[ind_celltype[:, np.newaxis], col_range] = ckdout[0]
    return knn_dists, knn_indices

def integrate(
          #input adata
          adata: AnnData, batch: str, cell_type: Optional[str] = None, use_rep: Optional[str] = None, n_latent: int = 50,
          #neighbors global setting
          n_neighbors: Optional[int] = None, n_meta_neighbors: int = 3, approx: bool = True, metric: Union[str, types.FunctionType, DistanceMetric] = 'euclidean',
          #if approx = True, annoy or pyNNDescent
          use_annoy: bool = True, annoy_n_trees: int = 10, pynndescent_n_neighbors: int = 30, pynndescent_random_state: int = 0,
          #if approx = False
          use_faiss: bool = True,
          #connectivities
          set_op_mix_ratio: float = 1.0, local_connectivity: int = 1, trim: Optional[int] = None,
          #random and copy
          neighbor_random_state: int = 0, copy: bool = False) -> Union[AnnData, None]:
    """
    Cell type controlled k nearest neighbors. This is a variant of BBKNN by searching neighbors across matched cell groups in different batches.
    For a given cell belonging to cell type 'c', first determine the batches that contain 'c' and its neighboring cell types, and then in each batch, search nearest neighbors out of them.

    Parameters
    ----------
    adata
        An :class:`~anndata.AnnData` object containing batch and cell type information in `.obs`, as well as latent space (e.g., `'X_pca'`) in `.obsm`.
    batch
        Column name (key) of cell metadata specifying batch information.
    cell_type
        Column name (key) of cell metadata specifying cell type information.
        Default to no cell type information provided (i.e., searching nearest neighbors in the entire batch space).
    use_rep
        Representation used to calculate distances. This can be any representations stored in `.obsm`.
        Default to the PCA coordinates (`'X_pca'`) if present.
    n_latent
        Number of latent representations used.
        Default to min(50, number of available latent representations).
    n_neighbors
        Total number of nearest neighbors for each cell. This number will be contributed equally from batches that qualify.
        Default to max(15, n) where n is the number of batches times three, meaning that each qualified batch will provide at least 3 neighbors.
        For example, if one cell type exists exclusively in one batch, then this batch needs to provide 15 neighbors.
    n_meta_neighbors
        Total number of nearest meta neighbors for each cell type in each batch (calculated from cell centroids).
        The final nearest meta neighbors are the union across batches that contain this given cell type.
        The smaller this value, the stronger bonding of the same cell type.
        Setting to 1 will make each cell search nearest neighbors only in the cell type it belongs to (i.e., forcibly clustering the same cell types).
        (Default: `3`)
    approx
        Whether to use fast approximate neighbor finding (annoy or pyNNDescent).
        (Default: `True`)
    metric
        Distance metric to use.
        (Default: `'euclidean'`)
    use_annoy
        Whether to use annoy for neighbor finding when `approx = True`. Setting `use_annoy = False` will use pyNNDescent instead.
        (Default: `True`)
    annoy_n_trees
        Number of trees to construct in the annoy forest when `approx = True` and `use_annoy = True`.
        (Default: `10`)
    pynndescent_n_neighbors
        Number of neighbors to include in the approximate neighbor graph when `approx = True` and `use_annoy = False`.
        (Default: `30`)
    pynndescent_random_state
        Random seed to use in pyNNDescent when `approx = True` and `use_annoy = False`.
        (Default: `0`)
    use_faiss
        Whether to use the faiss package to compute nearest neighbors if installed when `approx = False` and `metric = 'euclidean'`.
        (Default: `True`)
    set_op_mix_ratio
        Float between 0 and 1 controlling the blend between a connectivity matrix formed exclusively from mutual nearest neighbor pairs (0)
        and a union of all observed neighbor relationships with the mutual pairs emphasized (1).
        (Default: `1.0`)
    local_connectivity
        UMAP connectivity computation parameter controlling how many nearest neighbors of each cell are assumed to be fully connected (with a connectivity value of 1).
        (Default: `1`)
    trim
        Trim each cell to top `trim` connectivities. May help with population independence and improve the tidiness of clustering.
        Default to n_neighbors*10. Set to 0 to skip trimming.
    neighbor_random_state
        Random seed to use in assigning the remainder neighbors to batches.
        For example, assigning 10 nearest neighbors to 3 batches will make one remainder neighbor randomly assigned to one of the three batches.
        (Default: `0`)
    copy
        Whether to copy the adata or modify in-place.
        (Default: `False`)

    Returns
    ----------
    Union[AnnData, None]
        Depending on `copy`, return an updated or copied :class:`~anndata.AnnData` object with neighborhood graph included.
    """
    #check adata
    adata = adata.copy() if copy else adata
    if batch not in adata.obs:
        raise KeyError(
                f"🛑 '{batch}' is not found in the provided AnnData")
    batch_list = adata.obs[batch].astype(str).values
    if isinstance(cell_type, str) and cell_type not in adata.obs:
        raise KeyError(
                f"🛑 '{cell_type}' is not found in the provided AnnData")
    celltype_list = adata.obs[cell_type].astype(str).values if isinstance(cell_type, str) else np.full(adata.n_obs, 'cell', dtype = object)
    batch_counts = adata.obs[batch].astype(str).value_counts()
    few_batches = set(batch_counts.index[batch_counts <= 10])
    if len(few_batches) > 0:
        logger.warn(f"⚠️ The following batch(es) have too few cells (<= 10), please remove them before running `celltypist.integrate`: {few_batches}")
        return
    if use_rep is None:
        logger.info(f"👀 `use_rep` is not specified, will use `'X_pca'` as the search space")
        use_rep = 'X_pca'
    if use_rep not in adata.obsm.keys():
        raise KeyError(
                f"🛑 '{use_rep}' is not found in `.obsm`")
    n_latent = min([n_latent, adata.obsm[use_rep].shape[1]])
    pca = adata.obsm[use_rep][:, :n_latent]
    #check knn search params
    n_obs = adata.n_obs
    n_neighbors = max([15, 3 * len(np.unique(batch_list))]) if n_neighbors is None else n_neighbors
    swapped = False
    if approx:
        if use_annoy:
            computation = 'annoy'
            if metric not in ['angular', 'euclidean', 'manhattan', 'hamming']:
                swapped = True
                metric = 'euclidean'
        else:
            computation = 'pynndescent'
            if not (metric in pynndescent.distances.named_distances or isinstance(metric, types.FunctionType)):
                swapped = True
                metric = 'euclidean'
    else:
        if not ((metric == 'euclidean') or isinstance(metric, DistanceMetric) or metric in KDTree.valid_metrics):
            swapped = True
            metric = 'euclidean'
        if metric == 'euclidean':
            if 'faiss' in sys.modules and use_faiss:
                computation = 'faiss'
            else:
                computation = 'cKDTree'
        else:
            computation = 'KDTree'
    if swapped:
        logger.warn(f"👀 Unrecognized `metric` for type of neighbor calculation, will switch to 'euclidean'")
    #knn construction
    knn_dists, knn_indices = _get_graph(pca, batch_list, celltype_list, computation, n_neighbors, n_meta_neighbors, metric, annoy_n_trees, pynndescent_n_neighbors, pynndescent_random_state, neighbor_random_state)
    newidx = np.argsort(knn_dists, axis = 1)
    knn_indices = knn_indices[np.arange(n_obs)[:, np.newaxis], newidx]
    knn_dists = knn_dists[np.arange(n_obs)[:, np.newaxis], newidx]
    #connectivities + distances
    X = coo_matrix(([], ([], [])), shape=(n_obs, 1))
    connectivities = fuzzy_simplicial_set(X, n_neighbors, None, None, knn_indices = knn_indices, knn_dists = knn_dists, set_op_mix_ratio = set_op_mix_ratio, local_connectivity = local_connectivity)
    if isinstance(connectivities, tuple):
        connectivities = connectivities[0]
    connectivities = connectivities.tocsr()
    rows = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    cols = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    vals = np.zeros((n_obs * n_neighbors), dtype=np.float64)
    for i in range(n_obs):
        for j in range(n_neighbors):
            if knn_indices[i, j] == -1:
                continue
            if knn_indices[i, j] == i:
                val = 0.0
            else:
                val = knn_dists[i, j]
            rows[i * n_neighbors + j] = i
            cols[i * n_neighbors + j] = knn_indices[i, j]
            vals[i * n_neighbors + j] = val
    distances = coo_matrix((vals, (rows, cols)), shape=(n_obs, n_obs))
    distances.eliminate_zeros()
    #trim
    if trim is None:
        trim = 10 * n_neighbors
    if trim > 0:
        cutoffs = np.zeros(n_obs)
        for i in range(n_obs):
            row_array = connectivities.data[connectivities.indptr[i]: connectivities.indptr[i+1]]
            if row_array.shape[0] <= trim:
                continue
            cutoffs[i] = row_array[np.argsort(row_array)[-1*trim]]
        for iter in range(2):
            for i in range(n_obs):
                row_array = connectivities.data[connectivities.indptr[i]: connectivities.indptr[i+1]]
                row_array[row_array < cutoffs[i]] = 0
            connectivities.eliminate_zeros()
            connectivities = connectivities.T.tocsr()
    #assign
    adata.uns['neighbors'] = {}
    adata.uns['neighbors']['params'] = {'n_neighbors': n_neighbors, 'method': 'umap', 'metric': metric, 'n_pcs': n_latent, 'ccknn': {'n_meta_neighbors': n_meta_neighbors, 'trim': trim, 'computation': computation, 'batch': batch}}
    adata.uns['neighbors']['params']['use_rep'] = use_rep
    adata.obsp['distances'] = distances.tocsr()
    adata.obsp['connectivities'] = connectivities
    adata.uns['neighbors']['distances_key'] = 'distances'
    adata.uns['neighbors']['connectivities_key'] = 'connectivities'
    return adata if copy else None
