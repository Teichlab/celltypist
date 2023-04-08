from sklearn.metrics import pairwise_distances
import numpy as np
import pandas as pd
from scipy.sparse import spmatrix
from anndata import AnnData
from scipy.stats import rankdata
from typing import Optional
from .. import logger
from .symbols import NOVEL, REMAIN, UNASSIGN

class Distance():
    """
    Class that deals with the cross-dataset cell-by-cell-type distance matrix.

    Parameters
    ----------
    dist_mat
        Cell-by-cell-type distance matrix.
    cell
        Cell meta-information including at least `'dataset'`, `'ID'` and `'cell_type'`.
    cell_type
        Cell type meta-information including at least `'dataset'` and `'cell_type'`.

    Attributes
    ----------
    dist_mat
        A cell-by-cell-type distance matrix.
    cell
        Cell meta-information including `'dataset'`, `'ID'` and `'cell_type'`.
    cell_type
        Cell type meta-information including `'dataset'` and `'cell_type'`.
    n_cell
        Number of cells involved.
    n_cell_type
        Number of cell types involved.
    shape
        Tuple of number of cells and cell types.
    assignment
        Assignment of each cell to the most similar cell type in each dataset (obtained through the `assign` method).
    """
    def __init__(self, dist_mat: np.ndarray, cell: pd.DataFrame, cell_type: pd.DataFrame):
        self.dist_mat = dist_mat
        if cell.shape[0] != self.dist_mat.shape[0]:
            raise ValueError(
                    f"🛑 Number of cells in `cell` does not match the cell number in `dist_mat`")
        if cell_type.shape[0] != self.dist_mat.shape[1]:
            raise ValueError(
                    f"🛑 Number of cell types in `cell_type` does not match the cell type number in `dist_mat`")
        if not {'dataset', 'ID', 'cell_type'}.issubset(set(cell.columns)):
            raise KeyError(
                    f"🛑 Please include `'dataset'`, `'ID'` and `'cell_type'` as the cell meta-information")
        if not {'dataset', 'cell_type'}.issubset(set(cell_type.columns)):
            raise KeyError(
                    f"🛑 Please include `'dataset'` and `'cell_type'` as the cell type meta-information")
        self.cell = cell
        self.cell_type = cell_type

    @property
    def n_cell(self) -> int:
        """Number of cells."""
        return self.dist_mat.shape[0]

    @property
    def n_cell_type(self) -> int:
        """Number of cell types."""
        return self.dist_mat.shape[1]

    @property
    def shape(self) -> tuple:
        """Numbers of cells and cell types."""
        return self.dist_mat.shape

    def __repr__(self):
        lend = len(np.unique(self.cell_type.dataset))
        if lend > 1:
            base = f"Cross-dataset distance matrix between {self.n_cell} cells and {self.n_cell_type} cell types from {lend} datasets"
        else:
            base = f"Distance matrix between {self.n_cell} cells and {self.n_cell_type} cell types"
        base += f"\n    dist_mat: distance matrix between {self.n_cell} cells and {self.n_cell_type} cell types"
        base += f"\n    cell: cell meta-information ({str(list(self.cell.columns))[1:-1]})"
        base += f"\n    cell_type: cell type meta-information ({str(list(self.cell_type.columns))[1:-1]})"
        if hasattr(self, 'assignment'):
            base += f"\n    assignment: data frame of cross-dataset cell type assignment"
        return base

    @staticmethod
    def from_adata(adata: AnnData, dataset: str, cell_type: str, use_rep: Optional[str] = None, metric: Optional[str] = None, n_jobs: Optional[int] = None, check_params: bool = True, **kwargs):
        """
        Generate a :class:`~celltypist.contro.distance.Distance` object from the :class:`~anndata.AnnData` given.

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
            Number of CPUs used. Default to one CPU. `-1` means all CPUs are used.
        check_params
            Whether to check (or set the default) for `dataset`, `cell_type`, `use_rep` and `metric`.
            (Default: `True`)
        **kwargs
            Other keyword arguments passed to :func:`sklearn.metrics.pairwise_distances`.

        Returns
        ----------
        :class:`~celltypist.contro.distance.Distance`
            A :class:`~celltypist.contro.distance.Distance` object representing the cross-dataset cell-by-cell-type distance matrix.
        """
        #Use `check_params = False` if `dataset`, `cell_type`, `use_rep` and `metric` are already provided correctly.
        if check_params:
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
        Cell_X = adata.X if use_rep == 'X' else adata.obsm[use_rep]
        IDs = adata.obs_names
        datasets = adata.obs[dataset].astype(str).values
        celltypes = adata.obs[cell_type].astype(str).values
        use_Cell_X = Cell_X if use_rep != 'X' else np.expm1(Cell_X)
        Celltype_X = []
        col_ds = []
        col_cs =[]
        for d in np.unique(datasets):
            for c in np.unique(celltypes[datasets == d]):
                col_cs.append(c)
                col_ds.append(d)
                m = use_Cell_X[(datasets == d) & (celltypes == c), :].mean(axis = 0)
                Celltype_X.append(m.A1 if isinstance(m, np.matrix) else m)
        Celltype_X = np.log1p(np.array(Celltype_X)) if use_rep == 'X' else np.array(Celltype_X)
        if metric not in ['cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan']:
            if isinstance(Cell_X, spmatrix):
                Cell_X = Cell_X.toarray()
            if isinstance(Celltype_X, spmatrix):
                Celltype_X = Celltype_X.toarray()
        dist_mat = pairwise_distances(Cell_X, Celltype_X, metric = metric, n_jobs = n_jobs, **kwargs)
        cell = pd.DataFrame(dict(dataset=datasets, ID=IDs, cell_type=celltypes))
        cell_type = pd.DataFrame(dict(dataset=col_ds, cell_type=col_cs))
        return Distance(dist_mat, cell, cell_type)

    def normalize(self, Gaussian_kernel: bool = False, rank: bool = True, normalize: bool = True) -> None:
        """
        Normalize the distance matrix with a Gaussian kernel.

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
            The :class:`~celltypist.contro.distance.Distance` object modified with a normalized distance matrix.
        """
        if Gaussian_kernel:
            sds = np.sqrt((self.dist_mat ** 2).sum(axis = 1) / self.n_cell_type)[:, np.newaxis]
            self.dist_mat = np.exp(- self.dist_mat / (2 / sds)**2)
            self.dist_mat = 1 - self.dist_mat / self.dist_mat.sum(axis = 1)[:, np.newaxis]
        if rank:
            self.dist_mat = rankdata(self.dist_mat).reshape(self.dist_mat.shape)
        if normalize:
            self.dist_mat = self.dist_mat / self.dist_mat.max()

    def concatenate(self, *distances, by: str = 'cell', check: bool = False):
        """
        Concatenate by either cells (rows) or cell types (columns).

        Parameters
        ----------
        distances
            A :class:`~celltypist.contro.distance.Distance` object or a list of such objects.
        by
            The direction of concatenation, joining either cells (`'cell'`, rows) or cell types (`'cell_type'`, columns).
            (Default: `'cell'`)
        check
            Check whether the concatenation is feasible.
            (Default: `False`)

        Returns
        ----------
        :class:`~celltypist.contro.distance.Distance`
            A :class:`~celltypist.contro.distance.Distance` object concatenated along cells (`by = 'cell'`) or cell types (`by = 'cell_type'`).
        """
        distances = distances[0] if isinstance(distances[0], (list, tuple, set)) else distances
        distances = tuple(distances)
        all_distances = (self,) + distances
        if by not in ['cell', 'cell_type']:
            raise ValueError(
                    f"🛑 Unrecognized `by` value, should be one of `'cell'` or `'cell_type'`")
        if check:
            series_compare = [(x.cell_type.dataset+x.cell_type.cell_type).sort_values() for x in all_distances] if by == 'cell' else [(x.cell.dataset+x.cell.ID).sort_values() for x in all_distances]
            if pd.concat(series_compare, axis = 1).T.drop_duplicates().shape[0] > 1:
                raise Exception(
                        f"🛑 Concatenation is not feasible. Please ensure the meta-information is matched")
        if by == 'cell':
            dist_mat = np.concatenate([x.dist_mat for x in all_distances], axis = 0)
            cell = pd.concat([x.cell for x in all_distances], axis = 0, ignore_index = True)
            return Distance(dist_mat, cell, self.cell_type)
        else:
            match_base = (self.cell.dataset+self.cell.ID).reset_index().set_index(0)
            indices = [np.argsort(match_base.loc[x.cell.dataset+x.cell.ID, 'index'].values) for x in distances]
            dist_mat = np.concatenate([self.dist_mat] + [x.dist_mat[y, :] for x,y in zip(distances, indices)], axis = 1)
            cell_type = pd.concat([x.cell_type for x in all_distances], axis = 0, ignore_index = True)
            return Distance(dist_mat, self.cell, cell_type)

    def symmetric(self) -> bool:
        """
        Check whether the distance matrix is symmetric in terms of datasets and cell types.

        Returns
        ----------
        bool
            `True` or `False` indicating whether all datasets and cell types are included in the object (thus symmetric).
        """
        return np.array_equal(np.unique(self.cell.dataset + self.cell.cell_type), np.unique(self.cell_type.dataset + self.cell_type.cell_type))

    def filter_cells(self, check_symmetry: bool = True) -> None:
        """
        Filter out cells whose gene expression profiles do not correlate most with the eigen cell they belong to (i.e., correlate most with other cell types).

        Parameters
        ----------
        check_symmetry
            Whether to check the symmetry of the distance matrix in terms of datasets and cell types.
            (Default: `True`)

        Returns
        ----------
        None
            A :class:`~celltypist.contro.distance.Distance` object with undesirable cells filtered out.
        """
        if check_symmetry and not self.symmetric():
            raise ValueError(
                    f"🛑 Cell filtering is not possible. Please provide the matrix with symmetric datasets and cell types")
        bool_cell = np.ones(self.n_cell, dtype=bool)
        for i, s in self.cell.iterrows():
            flag_dataset = self.cell_type.dataset == s['dataset']
            if self.cell_type.cell_type.values[flag_dataset][self.dist_mat[i][flag_dataset].argmin()] != s['cell_type']:
                bool_cell[i] = False
        if (~bool_cell).sum() == 0:
            logger.info(f"✂️ No cells are filtered out")
        else:
            ds_unique, ds_table = np.unique(self.cell.dataset.values[~bool_cell], return_counts = True)
            if len(ds_unique) == 1:
                logger.info(f"✂️ {(~bool_cell).sum()} cells are filtered out from {ds_unique[0]}")
            else:
                logger.info(f"✂️ {(~bool_cell).sum()} cells are filtered out, including:")
                for m, n in zip(ds_unique, ds_table):
                    logger.info(f"      {n} cells from {m}")
            self.dist_mat = self.dist_mat[bool_cell]
            self.cell = self.cell[bool_cell]
            all_combine = (self.cell_type.dataset + ': ' + self.cell_type.cell_type).values
            left_combine = np.unique(self.cell.dataset + ': ' + self.cell.cell_type)
            if len(left_combine) < len(all_combine):
                column_keep = np.isin(all_combine, left_combine)
                self.dist_mat = self.dist_mat[:, column_keep]
                self.cell_type = self.cell_type[column_keep]
                logger.info(f"✂️ The following cell types are discarded due to low confidence in annotation:")
                for rec in all_combine[~column_keep]:
                    logger.info(f"      {rec}")

    def to_meta(self, check_symmetry: bool = True, turn_binary: bool = False, return_symmetry: bool = True) -> pd.DataFrame:
        """
        Meta-analysis of cross-dataset cell type dissimilarity or membership.

        Parameters
        ----------
        check_symmetry
            Whether to check the symmetry of the distance matrix in terms of datasets and cell types.
            (Default: `True`)
        turn_binary
            Whether to turn the distance matrix into a cell type membership matrix before meta analysis.
            (Default: `False`)
        return_symmetry
            Whether to return a symmetric dissimilarity matrix by averaging with its transposed form.
            (Default: `True`)

        Returns
        ----------
        :class:`~pandas.DataFrame`
            A :class:`~pandas.DataFrame` object representing the cell-type-level dissimilarity matrix (`turn_binary = False`) or membership matrix (`turn_binary = True`).
        """
        if check_symmetry and not self.symmetric():
            raise ValueError(
                    f"🛑 Meta cell analysis is not possible. Concatenate all datasets and cell types beforehand using `concatenate`")
        use_mat = self.to_binary(False if check_symmetry else True).dist_mat if turn_binary else self.dist_mat
        meta_cell = []
        for _, s in self.cell_type.iterrows():
            meta_cell.append(use_mat[(self.cell.dataset == s['dataset']) & (self.cell.cell_type == s['cell_type']), :].mean(axis = 0))
        meta_cell = pd.DataFrame(np.array(meta_cell))
        meta_cell.index = (self.cell_type.dataset + ': ' + self.cell_type.cell_type).values
        meta_cell.columns = meta_cell.index
        return (meta_cell + meta_cell.T)/2 if return_symmetry else meta_cell

    def to_binary(self, check_symmetry: bool = True):
        """
        Turn the distance matrix into a binary matrix representing the estimated cell type membership across datasets.

        Parameters
        ----------
        check_symmetry
            Whether to check the symmetry of the distance matrix in terms of datasets and cell types.
            (Default: `True`)

        Returns
        ----------
        :class:`~celltypist.contro.distance.Distance`
            A :class:`~celltypist.contro.distance.Distance` object representing the estimated cell type membership across datasets.
        """
        if check_symmetry and not self.symmetric():
            raise ValueError(
                    f"🛑 Cannot convert to a binary matrix. Please provide the matrix with symmetric datasets and cell types")
        member_mat = np.zeros(self.shape, dtype = int)
        datasets = self.cell_type.dataset.values
        for dataset in np.unique(datasets):
            indices = np.where(datasets == dataset)[0]
            member_mat[range(member_mat.shape[0]), indices[self.dist_mat[:, indices].argmin(axis = 1)]] = 1
        return Distance(member_mat, self.cell, self.cell_type)

    def assign(self) -> None:
        """
        Assign each cell to its most similar cell type in each dataset.

        Returns
        ----------
        None
            Modified object with the result of cell assignment added as `.assignment`.
        """
        assignment = {}
        for dataset in np.unique(self.cell_type.dataset):
            flag = self.cell_type.dataset == dataset
            assignment[dataset] = self.cell_type.cell_type.values[flag][self.dist_mat[:, flag].argmin(axis = 1)]
        assignment = pd.DataFrame(assignment, index = self.cell.index)
        #no need to assign cells for the dataset they belong to
        for dataset in assignment.columns:
            flag = self.cell.dataset == dataset
            assignment.loc[flag, dataset] = self.cell.cell_type.values[flag]
        self.assignment = assignment

    def to_confusion(self, D1: str, D2: str, check: bool = True) -> tuple:
        """
        This function is deprecated. Use `to_pairwise_confusion` and `to_multi_confusion` instead.
        Extract the dataset1-by-dataset2 and dataset2-by-dataset1 confusion matrices. Note this function is expected to be applied to a binary membership matrix.

        Parameters
        ----------
        D1
            Name of the first dataset.
        D2
            Name of the second dataset.
        check
            Whether to check names of the two datasets are contained.
            (Default: `True`)

        Returns
        ----------
        tuple
            The dataset1-by-dataset2 and dataset2-by-dataset1 confusion matrices.
        """
        if check and not {D1, D2}.issubset(np.unique(self.cell_type.dataset)):
            raise ValueError(
                    f"🛑 Please provide correct dataset names")
        D1_col_flag = self.cell_type.dataset == D1
        D2_col_flag = self.cell_type.dataset == D2
        D1_celltypes = self.cell_type.cell_type.values[D1_col_flag]
        D2_celltypes = self.cell_type.cell_type.values[D2_col_flag]
        D1_row_flag = self.cell.dataset == D1
        D2_row_flag = self.cell.dataset == D2
        D1byD2 = pd.DataFrame(np.array([self.dist_mat[D1_row_flag & (self.cell.cell_type == x)][:, D2_col_flag].sum(axis=0) for x in D1_celltypes]), columns = D2_celltypes, index = D1_celltypes)
        D2byD1 = pd.DataFrame(np.array([self.dist_mat[D2_row_flag & (self.cell.cell_type == x)][:, D1_col_flag].sum(axis=0) for x in D2_celltypes]), columns = D1_celltypes, index = D2_celltypes)
        return D1byD2, D2byD1

    def to_pairwise_confusion(self, D1: str, D2: str, check: bool = True) -> tuple:
        """
        Extract the dataset1-by-dataset2 and dataset2-by-dataset1 confusion matrices.

        Parameters
        ----------
        D1
            Name of the first dataset.
        D2
            Name of the second dataset.
        check
            Whether to check names of the two datasets are contained.
            (Default: `True`)

        Returns
        ----------
        tuple
            The dataset1-by-dataset2 and dataset2-by-dataset1 confusion matrices.
        """
        if check and not {D1, D2}.issubset(np.unique(self.cell_type.dataset)):
            raise ValueError(
                    f"🛑 Please provide correct dataset names")
        if not hasattr(self, 'assignment'):
            raise AttributeError(
                    f"🛑 No `.assignment` attribute in the object. Use the `.assign` method first")
        D1_flag = (self.cell.dataset == D1)
        D2_flag = (self.cell.dataset == D2)
        D1byD2 = pd.crosstab(self.cell.cell_type[D1_flag], self.assignment.loc[D1_flag, D2])
        D2byD1 = pd.crosstab(self.cell.cell_type[D2_flag], self.assignment.loc[D2_flag, D1])
        D1byD2_lack_columns = D2byD1.index.difference(D1byD2.columns)
        if len(D1byD2_lack_columns) > 0:
            D1byD2 = D1byD2.join(pd.DataFrame(np.zeros((len(D1byD2.index), len(D1byD2_lack_columns)), dtype=int), index = D1byD2.index, columns = D1byD2_lack_columns))
        D2byD1_lack_columns = D1byD2.index.difference(D2byD1.columns)
        if len(D2byD1_lack_columns) > 0:
            D2byD1 = D2byD1.join(pd.DataFrame(np.zeros((len(D2byD1.index), len(D2byD1_lack_columns)), dtype=int), index = D2byD1.index, columns = D2byD1_lack_columns))
        return D1byD2, D2byD1.loc[D1byD2.columns, D1byD2.index]

    def to_multi_confusion(self, relation: pd.DataFrame, D: str, check: bool = True) -> tuple:
        """
        Extract the confusion matrices between meta-cell-types defined prior and cell types from a new dataset.

        Parameters
        ----------
        relation
            A :class:`~pandas.DataFrame` object representing the cell type harmonization result across multiple datasets.
        D
            Name of the new dataset to be aligned.
        check
            Whether to check names of the datasets are contained.
            (Default: `True`)

        Returns
        ----------
        tuple
            The confusion matrices between meta-cell-types defined prior and cell types from a new dataset.
        """
        datasets = relation.columns[0::2]
        if check:
            if not set(datasets).issubset(np.unique(self.cell_type.dataset)):
                raise ValueError(
                        f"🛑 `relation` contains unexpected dataset names")
            if D not in np.unique(self.cell_type.dataset) or D in datasets:
                raise ValueError(
                        f"🛑 Please provide a valid dataset name `D`")
        if not hasattr(self, 'assignment'):
            raise AttributeError(
                    f"🛑 No `.assignment` attribute in the object. Use the `.assign` method first")
        #D1byD2
        D1_flag = self.cell.dataset.isin(datasets)
        D1_assign = self.assignment[D1_flag]
        D1_truth = np.full(D1_assign.shape[0], UNASSIGN, dtype = object)
        for _, s in relation.iterrows():
            celltypes = s.values[0::2]
            non_blank_flag = ~np.isin(celltypes, [NOVEL, REMAIN])
            existing_datasets = datasets[non_blank_flag]
            existing_celltypes = celltypes[non_blank_flag]
            flag = np.all(D1_assign[existing_datasets] == existing_celltypes, axis = 1).values & self.cell[D1_flag].dataset.isin(existing_datasets).values
            D1_truth[flag] = ' '.join(s.values)
        D1_used = D1_truth != UNASSIGN
        D1byD2 = pd.crosstab(D1_truth[D1_used], D1_assign.loc[D1_used, D])
        #D2byD1
        D2_flag = self.cell.dataset == D
        D2_assign = self.assignment[D2_flag]
        D2_predict = np.full(D2_assign.shape[0], UNASSIGN, dtype = object)
        for _, s in relation.iterrows():
            celltypes = s.values[0::2]
            flags = (D2_assign[datasets] == celltypes) | np.isin(celltypes, [NOVEL, REMAIN])
            D2_predict[np.all(flags, axis = 1).values] = ' '.join(s.values)
        D2_used = D2_predict != UNASSIGN
        D2byD1 = pd.crosstab(self.cell.cell_type[D2_flag][D2_used], D2_predict[D2_used])
        #warning
        if relation.shape[0] > D1byD2.shape[0]:
            lost_celltypes = np.setdiff1d(relation.apply(lambda row: ' '.join(row.values), axis = 1).values, D1byD2.index)
            logger.warn(f"⚠️ Warning: no cells are found to match these patterns: {set(lost_celltypes)}. Double check the harmonized relationships before integrating '{D}'")
            D1byD2 = pd.concat([D1byD2, pd.DataFrame(np.zeros((len(lost_celltypes), len(D1byD2.columns)), dtype=int), index = lost_celltypes, columns = D1byD2.columns)], axis = 0)
        #a unique cell type in D2 may be annotated to nothing and filtered
        lost_celltypes = np.setdiff1d(np.unique(self.cell.cell_type[D2_flag]), D2byD1.index)
        if len(lost_celltypes) > 0:
            D2byD1 = pd.concat([D2byD1, pd.DataFrame(np.zeros((len(lost_celltypes), len(D2byD1.columns)), dtype=int), index = lost_celltypes, columns = D2byD1.columns)], axis = 0)
        #return
        D1byD2_lack_columns = D2byD1.index.difference(D1byD2.columns)
        if len(D1byD2_lack_columns) > 0:
            D1byD2 = D1byD2.join(pd.DataFrame(np.zeros((len(D1byD2.index), len(D1byD2_lack_columns)), dtype=int), index = D1byD2.index, columns = D1byD2_lack_columns))
        D2byD1_lack_columns = D1byD2.index.difference(D2byD1.columns)
        if len(D2byD1_lack_columns) > 0:
            D2byD1 = D2byD1.join(pd.DataFrame(np.zeros((len(D2byD1.index), len(D2byD1_lack_columns)), dtype=int), index = D2byD1.index, columns = D2byD1_lack_columns))
        return D1byD2, D2byD1.loc[D1byD2.columns, D1byD2.index]
