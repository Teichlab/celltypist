import os
from anndata import AnnData
from typing import Optional, Union
import numpy as np
import pandas as pd

_samples_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "samples")


def _get_sample_data(filename: str) -> str:
    """Get the full path to the sample input data included in the package."""
    return os.path.join(_samples_path, filename)


def get_sample_csv() -> str:
    """
    Get the full path to the sample csv file included in the package.

    Returns
    ----------
    str
        A string of the full path to the sample csv file (`sample_cell_by_gene.csv`).
    """
    return _get_sample_data("sample_cell_by_gene.csv")

def downsample_adata(adata: AnnData,
                     mode: str = 'total',
                     n_cells: Optional[int] = None,
                     by: Optional[str] = None,
                     balance_cell_type: bool = False,
                     random_state: int = 0,
                     return_index: bool = True) -> Union[AnnData, np.ndarray]:
    """
    Downsample cells to a given number (either in total or per cell type).

    Parameters
    ----------
    adata
        An :class:`~anndata.AnnData` object representing the input data.
    mode
        The way downsampling is performed. Default to downsampling the input cells to a total of `n_cells`.
        Set to `'each'` if you want to downsample cells within each cell type to `n_cells`.
        (Default: `'total'`)
    n_cells
        The total number of cells (`mode = 'total'`) or the number of cells from each cell type (`mode = 'each'`) to sample.
        For the latter, all cells from a given cell type will be selected if its cell number is fewer than `n_cells`.
    by
        Key (column name) of the input AnnData representing the cell types.
    balance_cell_type
        Whether to balance the cell type frequencies when `mode = 'total'`.
        Setting to `True` will sample rare cell types with a higher probability, ensuring close-to-even cell type compositions.
        This argument is ignored if `mode = 'each'`.
        (Default: `False`)
    random_state
        Random seed for reproducibility.
    return_index
        Only return the downsampled cell indices.
        Setting to `False` if you want to get a downsampled version of the input AnnData.
        (Default: `True`)

    Returns
    ----------
    Depending on `return_index`, returns the downsampled cell indices or a subset of the input AnnData.
    """
    np.random.seed(random_state)
    if n_cells is None:
        raise ValueError(
                f"🛑 Please provide `n_cells`")
    if mode == 'total':
        if n_cells >= adata.n_obs:
            raise ValueError(
                    f"🛑 `n_cells` ({n_cells}) should be fewer than the total number of cells ({adata.n_obs})")
        if balance_cell_type:
            if by is None:
                raise KeyError(
                        f"🛑 Please specify the cell type column if you want to balance the cell type frequencies")
            labels = adata.obs[by]
            celltype_freq = np.unique(labels, return_counts = True)
            len_celltype = len(celltype_freq[0])
            mapping = pd.Series(1 / (celltype_freq[1]*len_celltype), index = celltype_freq[0])
            p = mapping[labels].values
            sampled_cell_index = np.random.choice(adata.n_obs, n_cells, replace = False, p = p)
        else:
            sampled_cell_index = np.random.choice(adata.n_obs, n_cells, replace = False)
    elif mode == 'each':
        if by is None:
            raise KeyError(
                    f"🛑 Please specify the cell type column for downsampling")
        celltypes = np.unique(adata.obs[by])
        sampled_cell_index = np.concatenate([np.random.choice(np.where(adata.obs[by] == celltype)[0], min([n_cells, np.sum(adata.obs[by] == celltype)]), replace = False) for celltype in celltypes])
    else:
        raise ValueError(
                f"🛑 Unrecognized `mode` value, should be one of `'total'` or `'each'`")
    if return_index:
        return sampled_cell_index
    else:
        return adata[sampled_cell_index].copy()
