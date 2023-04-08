from anndata import AnnData
from typing import Union, Optional
import pandas as pd
import numpy as np
from .distance import Distance
from .pct import PredictiveClusteringTree
from .. import logger

class Distances():
    """
    Class that deals with the cross-dataset cell and cell type distance/dissimilarity analysis.
    For each dataset, a cell-by-cell-type distance matrix is first calculated. This matrix is propagated to other datasets (with the same set of cell types) using predictive clustering trees. Finally, matrices are concatenated across all datasets.

    Parameters
    ----------
    adata
        An :class:`~anndata.AnnData` object containing different datasets/batches and cell types.
        In most scenarios, the format of the expression `.X` in the AnnData is flexible (normalized, log-normalized, z-scaled, etc.).
        However, when `use_rep` is specified as `'X'` (or `X_pca` is not detected in `.obsm` and no other latent representations are provided), `.X` should be log-normalized (to a constant total count per cell).
    dataset
        Column name (key) of cell metadata specifying dataset information.
    cell_type
        Column name (key) of cell metadata specifying cell type information.
    use_rep
        Representation used to calculate distances. This can be `'X'` or any representations stored in `.obsm`.
        Default to the PCA coordinates if present (if not, use the expression matrix `X`).
    metric
        Metric to calculate the distance between each cell and each cell type. Can be `'euclidean'`, `'cosine'`, `'manhattan'` or any metrics applicable to :func:`sklearn.metrics.pairwise_distances`.
        Default to `'euclidean'` if latent representations are used for calculating distances, and to `'correlation'` if the expression matrix is used.
    n_jobs
        Number of CPUs used to calculate intra-dataset distances between cells and cell types. Default to one CPU. `-1` means all CPUs are used.
    **kwargs
        Other keyword arguments passed to :func:`sklearn.metrics.pairwise_distances`.

    Attributes
    ----------
    datasets
        List of datasets involved.
    distances
        List of :class:`~celltypist.contro.distance.Distance` objects, with the ith element for dataset i.
    cell_types
        List of cell types. The ith element lists the cell types in dataset i.
    cells
        List of cells. The ith element lists cells in dataset i.
    adata
        The AnnData that is referred.
    pcts
        List of :class:`~celltypist.contro.pct.PredictiveClusteringTree` objects (after training), with the ith element for dataset i.
    """
    def __init__(self, adata: AnnData, dataset: str, cell_type: str, use_rep: Optional[str] = None, metric: Optional[str] = None, n_jobs: Optional[int] = None, **kwargs):
        if dataset not in adata.obs:
            raise KeyError(
                    f"🛑 '{dataset}' is not found in the provided AnnData")
        if cell_type not in adata.obs:
            raise KeyError(
                    f"🛑 '{cell_type}' is not found in the provided AnnData")
        if use_rep is None:
            if 'X_pca' in adata.obsm.keys():
                logger.info(f"👀 Detected PCA coordinates in the object, will use these to calculate distances")
                use_rep = 'X_pca'
            else:
                logger.info(f"🧙 Using the expression matrix to calculate distances")
                use_rep = 'X'
        elif (use_rep not in adata.obsm.keys()) and (use_rep != 'X'):
            raise KeyError(
                    f"🛑 '{use_rep}' is not found in `.obsm`")
        if use_rep == 'X' and adata.n_vars > 15000:
            logger.warn(f"⚠️ Warning: {adata.n_vars} features are used for calculating distances. Subsetting the AnnData into HVGs is recommended")
        if metric is None:
            metric = 'correlation' if use_rep == 'X' else 'euclidean'
        self.datasets = np.unique(adata.obs[dataset])
        self.distances = []
        for d in self.datasets:
            self.distances.append(Distance.from_adata(adata[adata.obs[dataset] == d], dataset, cell_type, use_rep, metric, n_jobs, False, **kwargs))
        self.use_rep = use_rep
        self.metric = metric
        self.adata = adata

    @property
    def cell_types(self) -> list:
        """Get the cell type list."""
        return [x.cell_type.cell_type.values for x in self.distances]

    @property
    def cells(self) -> list:
        """Get the cell list."""
        return [x.cell.ID.values for x in self.distances]

    def __repr__(self):
        base = f"Within-dataset distance matrices for {len(self.datasets)} datasets"
        base += f"\n    datasets: {str(list(self.datasets))[1:-1]}"
        base += f"\n    distances: list of distances for each dataset"
        if hasattr(self, 'pcts'):
            base += f"\n    pcts: list of predictive clustering trees for each dataset"
        base += f"\n    adata: AnnData object referred"
        return base

    def filter_cells(self, check_symmetry: bool = False) -> None:
        """
        For each dataset, filter out cells whose gene expression profiles do not correlate most with the eigen cell they belong to (i.e., correlate most with other cell types).

        Parameters
        ----------
        check_symmetry
            Whether to check the symmetry of the distance matrix in terms of datasets and cell types.
            (Default: `False`)

        Returns
        ----------
        None
           Modified :class:`~celltypist.contro.distances.Distances` object with undesirable cells filtered out.
        """
        for i in range(len(self.datasets)):
            self.distances[i].filter_cells(check_symmetry)

    def normalize(self, Gaussian_kernel: bool = False, rank: bool = True, normalize: bool = True) -> None:
        """
        Normalize the distance matrix for each dataset with a Gaussian kernel.

        Parameters
        ----------
        Gaussian_kernel
            Whether to apply the Gaussian kernel to the distance matrix.
            (Default: `False`)
        rank
            Whether to turn the matrix into a rank matrx.
            (Default: `True`)
        normalize
            Whether to maximum-normalize the distance matrix.
            (Default: `True`)

        Returns
        ----------
        None
            Modified :class:`~celltypist.contro.distances.Distances` object with normalized intra-dataset distance matrices.
        """
        for i in range(len(self.datasets)):
            self.distances[i].normalize(Gaussian_kernel = Gaussian_kernel, rank = rank, normalize = normalize)

    def train(self,
                max_depth: Optional[int] = None,
                min_samples_split: Union[int, float] = 20,
                min_samples_leaf: Union[int, float] = 10,
                min_weight_fraction_leaf: float = 0.0,
                random_state: Optional[int] = None,
                max_leaf_nodes: Optional[int] = None,
                F_test_prune: bool = True,
                p_thres: float = 0.05,
                sample_weight = None) -> None:
        """
        Train a predictive clustering tree (PCT) for each dataset using the intra-dataset distance matrix.

        Parameters
        ----------
        max_depth
            Maximum possible depth of the tree, starting from the root node which has a depth of 0.
            Default to no limit.
        min_samples_split
            The minimum sample size (in absolute number or fraction) of a possible node.
            (Default: `20`)
        min_samples_leaf
            The minimum sample size (in absolute number or fraction) of a possible leaf.
            (Default: `10`)
        min_weight_fraction_leaf
            The minimum fraction out of total sample weights for a possible leaf.
            (Default: `0.0`)
        random_state
            Random seed for column (feature) shuffling before selecting the best feature and threshold.
        max_leaf_nodes
            The maximum number of leaves, achieved by keeping high-quality (i.e., high impurity reduction) nodes.
            Default to no limit.
        F_test_prune
            Whether to use a F-test to prune the tree by removing unnecessary splits.
            (Default: `True`)
        p_thres
            p-value threshold for pruning nodes after F-test.
            (Default: `0.05`)
        sample_weight
            Sample weights. Default to equal weights across cells.

        Returns
        ----------
        None
            Modified :class:`~celltypist.contro.distances.Distances` object with PCT regressors added as `.pcts`.
        """
        self.pcts = []
        logger.info(f"🏋️ Training the predictive clustering trees for:")
        for da,c,di in zip(self.datasets, self.cells, self.distances):
            logger.info(f"      {da}")
            PCT = PredictiveClusteringTree(max_depth = max_depth,
                                           min_samples_split = min_samples_split,
                                           min_samples_leaf = min_samples_leaf,
                                           min_weight_fraction_leaf = min_weight_fraction_leaf,
                                           random_state = random_state,
                                           max_leaf_nodes = max_leaf_nodes,
                                           F_test_prune = F_test_prune,
                                           p_thres = p_thres)
            PCT.fit(self.adata[c].X, di.dist_mat, sample_weight = sample_weight)
            self.pcts.append(PCT)

    def predict(self, normalize: bool = True, return_distance: bool = True, **kwargs) -> Union[Distance, list]:
        """
        Iteratively predict the cell-by-cell-type distance matrix in other datasets using the predictive clustering tree (PCT) from one dataset.

        Parameters
        ----------
        normalize
            Whether to normalize the distance for each dataset.
            (Default: `True`)
        return_distance
            Whether to join the distances across datasets after prediction to form an entire :class:`~celltypist.contro.distance.Distance` object.
            (Default: `True`)
        **kwargs
            Keyword arguments passed to :func:`~celltypist.contro.distance.Distance.normalize`.

        Returns
        ----------
        If `return_distance` is `True`, a :class:`~celltypist.contro.distance.Distance` object will be returned representing all datasets and cell types combined; otherwise, a list will be returned, with the ith element for dataset i.
        """
        logger.info(f"🖋️ Predicting distances using the PCT in each dataset")
        new_ds = []
        for i in range(len(self.datasets)):
            dist_mat = self.pcts[i].predict(self.adata[np.concatenate(self.cells[:i] + self.cells[i+1:])].X)
            cell = pd.concat([x.cell for x in self.distances[:i]+self.distances[i+1:]], axis = 0)
            pred_dis = Distance(dist_mat, cell, self.distances[i].cell_type)
            new_ds.append(self.distances[i].concatenate(pred_dis, by = 'cell'))
        if normalize:
            for each_new_ds in new_ds:
                each_new_ds.normalize(**kwargs)
        if return_distance:
            return new_ds[0].concatenate(new_ds[1:], by = 'cell_type')
        else:
            return new_ds
