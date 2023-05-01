from anndata import AnnData
from typing import Union, Optional
import pandas as pd
import numpy as np
from .distance import Distance
from .distances import Distances
from .align import DistanceAlignment
from .. import logger

def harmonize(adata: AnnData,
              #dataset & cell type info
              dataset: str, cell_type: str, use_rep: Optional[str] = None, metric: Optional[str] = None,
              #PCT the distances or directly calculate
              use_pct: bool = False,
              #filter and normalize; note PCT will normalize after prediction
              filter_cells: bool = False, normalize: bool = True, Gaussian_kernel: bool = False,
              #PCT train
              F_test_prune: bool = True, p_thres: float = 0.05, random_state: int = 2,
              #order of datasets
              dataset_order: Optional[Union[list, tuple, np.ndarray, pd.Series, pd.Index]] = None, reorder_dataset: bool = True,
              #align params
              minimum_unique_percents: Union[list, tuple, np.ndarray, pd.Series, pd.Index, float] = (0.4, 0.5, 0.6, 0.7, 0.8),
              minimum_divide_percents: Union[list, tuple, np.ndarray, pd.Series, pd.Index, float] = (0.1, 0.15, 0.2),
              maximum_novel_percent: float = 0.05,
              #reannotate
              reannotate: bool = True, add_group: bool = True, prefix: str = '',
              #to PCT train
              **kwargs) -> DistanceAlignment:
    """
    PCT-based cell type harmonization across datasets/batches.

    Parameters
    ----------
    adata
        An :class:`~anndata.AnnData` object containing different datasets/batches and cell types.
        In most scenarios, the format of the expression `.X` in the AnnData is flexible (normalized, log-normalized, z-scaled, etc.).
        However, when `use_rep` is specified as `'X'` (or `X_pca` is not detected in `.obsm` and no other latent representations are provided), `.X` should be log-normalized (to a constant total count per cell).
    dataset
        Column name (key) of cell metadata specifying dataset/batch information.
    cell_type
        Column name (key) of cell metadata specifying cell type information.
    use_rep
        Representation used to calculate distances. This can be `'X'` or any representations stored in `.obsm`.
        Default to the PCA coordinates if present (if not, use the expression matrix `X`).
    metric
        Metric to calculate the distance between each cell and each cell type. Can be `'euclidean'`, `'cosine'`, `'manhattan'` or any metrics applicable to :func:`sklearn.metrics.pairwise_distances`.
        Default to `'euclidean'` if latent representations are used for calculating distances, and to `'correlation'` if the expression matrix is used.
    use_pct
        Whether to use a predictive clustering tree to infer cross-dataset cell type distances.
        Setting to `True` will calculate distances based on PCT, which is intended for datasets with large batch effects.
        (Default: `False`)
    filter_cells
        Whether to filter out cells whose gene expression profiles do not correlate most with the eigen cell they belong to (i.e., correlate most with other cell types).
        Setting to `True` will speed up the run as only a subset of cells are used, but will render the remaining cells (i.e., filtered cells) unannotated (see the `reannotate` argument).
        (Default: `False`)
    normalize
        Whether to normalize the distance matrix if `use_pct = False` (or normalize the predicted distance if `use_pct = True`).
        (Default: `True`)
    Gaussian_kernel
        Whether to apply the Gaussian kernel to the distance matrix.
        (Default: `False`)
    F_test_prune
        Whether to use a F-test to prune the tree by removing unnecessary splits.
        (Default: `True`)
    p_thres
        p-value threshold for pruning nodes after F-test.
        (Default: `0.05`)
    random_state
        Random seed for feature shuffling during PCT training.
        (Default: `2`)
    dataset_order
        Order of datasets to be aligned. If this argument is specified, `reorder_dataset` is ignored.
        Default to the order in the distance matrix (alphabetical order in most cases) if `reorder_dataset = False`.
    reorder_dataset
        Whether to reorder datasets based on their pairwise similarities.
        (Default: `True`)
    minimum_unique_percents
        The minimum cell assignment fraction(s) to claim a cell type as uniquely matched to a cell type from the other dataset.
        By default, five values will be tried (0.4, 0.5, 0.6, 0.7, 0.8) to find the one that produces least alignments in each harmonization iteration.
    minimum_divide_percents
        The minimum cell assignment fraction(s) to claim a cell type as divisible into two or more cell types from the other dataset.
        By default, three values will be tried (0.1, 0.15, 0.2) to find the one that produces least alignments in each harmonization iteration.
    maximum_novel_percent
        The maximum cell assignment fraction to claim a cell type as novel to a given dataset.
        (Default: `0.05`)
    reannotate
        Whether to reannotate cells into harmonized cell types.
        (Default: `True`)
    add_group
        Whether to annotate out cell type group information as well during reannotation.
        (Default: `True`)
    prefix
        Column prefix for the reannotation data frame.
    **kwargs
        Other keyword arguments passed to :class:`~celltypist.contro.pct.PredictiveClusteringTree`.

    Returns
    ----------
    DistanceAlignment
        A :class:`~celltypist.contro.align.DistanceAlignment` object. Four important attributes within this class are:
        1) :attr:`~celltypist.contro.align.DistanceAlignment.base_distance`, cross-dataset distances between all cells and all cell types.
        2) :attr:`~celltypist.contro.align.DistanceAlignment.relation`, the harmonization table.
        3) :attr:`~celltypist.contro.align.DistanceAlignment.groups`, high-hierarchy cell types categorizing rows of the harmonization table.
        4) :attr:`~celltypist.contro.align.DistanceAlignment.reannotation`, reannotated cell types and cell type groups.
    """
    #raw counts are not allowed to build trees
    if use_pct and adata.X.min() >= 0 and float(adata.X.max()).is_integer():
        raise ValueError(
                f"🛑 `.X` of the AnnData is detected to be raw counts, which is not suitable for building PCT")
    #build PCT using all genes is not realistic
    if use_pct and adata.n_vars > 15000:
        logger.warn(f"⚠️ Warning: {adata.n_vars} features are used and may take long time for building PCT. Subsetting the AnnData into HVGs is recommended")
    #generate a combined `Distance`
    if use_pct:
        separate_distances = Distances(adata, dataset = dataset, cell_type = cell_type, use_rep = use_rep, metric = metric, n_jobs = -1)
        if filter_cells:
            separate_distances.filter_cells(check_symmetry = False)
        separate_distances.train(F_test_prune = F_test_prune, p_thres = p_thres, random_state = random_state, **kwargs)
        combined_distance = separate_distances.predict(normalize = normalize, return_distance = True, Gaussian_kernel = Gaussian_kernel)
    else:
        combined_distance = Distance.from_adata(adata, dataset = dataset, cell_type = cell_type, use_rep = use_rep, metric = metric, n_jobs = -1, check_params = True)
        if filter_cells:
            combined_distance.filter_cells(check_symmetry = False)
        if normalize:
            combined_distance.normalize(Gaussian_kernel = Gaussian_kernel, rank = True, normalize = True)
    #before cell type alignment
    combined_distance.assign()
    alignment = DistanceAlignment(combined_distance, check = False, dataset_order = dataset_order, row_normalize = True, maximum_novel_percent = maximum_novel_percent)
    if dataset_order is None and reorder_dataset:
        logger.info(f"🏆 Reordering datasets")
        alignment.reorder_dataset()
    #cell type alignment
    alignment.best_align(dataset_order = None, minimum_unique_percents = minimum_unique_percents, minimum_divide_percents = minimum_divide_percents)
    #reannotate
    if reannotate:
        logger.info(f"🖋️ Reannotating cells")
        alignment.reannotate(show_iteration = False, add_group = add_group, prefix = prefix)
    logger.info(f"✅ Harmonization done!")
    #return
    return alignment

harmonise = harmonize
