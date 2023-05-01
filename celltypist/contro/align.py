import numpy as np
import pandas as pd
from typing import Union, Optional
from .. import logger
from .distance import Distance
from .symbols import ONE2ONE, ONE2MANY, MANY2ONE, NOVEL, REMAIN, UNASSIGN
from .plot import SEP1, _identify_relation_groups
import os
import pickle

def _subset_and_normalize(D1byD2: pd.DataFrame, D2byD1: pd.DataFrame, relation: pd.DataFrame, row_normalize: bool = True) -> tuple:
    """
    For internal use. Subset the two matrices and row normalize to a sum of 1.
    """
    if relation.shape[0] > 0:
        flag1 = ~D1byD2.index.isin(relation.D1.unique())
        flag2 = ~D1byD2.columns.isin(relation.D2.unique())
        D1byD2 = D1byD2.loc[flag1, flag2]
        D2byD1 = D2byD1.loc[flag2, flag1]
    if not D1byD2.empty and row_normalize:
        D1byD2 = D1byD2 / D1byD2.sum(axis = 1).values[:, np.newaxis]
        D2byD1 = D2byD1 / D2byD1.sum(axis = 1).values[:, np.newaxis]
    return D1byD2, D2byD1

def _pairwise_align(D1byD2: pd.DataFrame, D2byD1: pd.DataFrame, check: bool = True, row_normalize: bool = True, minimum_unique_percent: float = 0.5, minimum_divide_percent: float = 0.1, maximum_novel_percent: float = 0.05) -> pd.DataFrame:
    """
    For internal use. Align cell types from two datasets based on two confusion matrices.
    """
    if check and (not np.array_equal(D1byD2.index, D2byD1.columns) or not np.array_equal(D1byD2.columns, D2byD1.index)):
        raise ValueError(
                f"🛑 Matrices are not provided correctly")
    relation = pd.DataFrame(columns = ['D1', 'relation', 'D2'])
    #rows all 0
    rownames = D1byD2.index[D1byD2.sum(axis = 1) == 0]
    if len(rownames) > 0:
        relation = pd.concat([relation, pd.DataFrame({'D1': rownames, 'relation': [ONE2ONE]*len(rownames), 'D2': [NOVEL]*len(rownames)})])
    rownames = D2byD1.index[D2byD1.sum(axis = 1) == 0]
    if len(rownames) > 0:
        relation = pd.concat([relation, pd.DataFrame({'D1': [NOVEL]*len(rownames), 'relation': [ONE2ONE]*len(rownames), 'D2': rownames})])
    #subset and normalize
    D1byD2, D2byD1 = _subset_and_normalize(D1byD2, D2byD1, relation, row_normalize)
    #D1 novel
    colnames = D2byD1.columns[D2byD1.max(axis = 0) < maximum_novel_percent]
    if len(colnames) > 0:
        relation = pd.concat([relation, pd.DataFrame({'D1': colnames, 'relation': [ONE2ONE]*len(colnames), 'D2': [NOVEL]*len(colnames)})])
    #D2 novel
    colnames = D1byD2.columns[D1byD2.max(axis = 0) < maximum_novel_percent]
    if len(colnames) > 0:
        relation = pd.concat([relation, pd.DataFrame({'D1': [NOVEL]*len(colnames), 'relation': [ONE2ONE]*len(colnames), 'D2': colnames})])
    D1byD2, D2byD1 = _subset_and_normalize(D1byD2, D2byD1, relation, row_normalize)
    #one2one
    indices = (D1byD2 > minimum_unique_percent) & (D2byD1.T > minimum_unique_percent)
    one2many_indices = (D1byD2 > minimum_divide_percent) & (D2byD1.T > minimum_unique_percent)
    one2many_indices = one2many_indices & (one2many_indices.sum(axis = 1) >= 2).values[:, np.newaxis]
    many2one_indices = (D2byD1 > minimum_divide_percent) & (D1byD2.T > minimum_unique_percent)
    many2one_indices = many2one_indices & (many2one_indices.sum(axis = 1) >= 2).values[:, np.newaxis]
    row_index, col_index = np.where(indices & (indices.sum(axis = 1) == 1).values[:, np.newaxis] & (indices.sum(axis = 0) == 1).values & (~one2many_indices) & (~many2one_indices.T))
    while len(row_index) > 0:
        relation = pd.concat([relation, pd.DataFrame({'D1': D1byD2.index[row_index], 'relation': [ONE2ONE]*len(row_index), 'D2': D1byD2.columns[col_index]})])
        D1byD2, D2byD1 = _subset_and_normalize(D1byD2, D2byD1, relation, row_normalize)
        if D1byD2.empty:
            break
        indices = (D1byD2 > minimum_unique_percent) & (D2byD1.T > minimum_unique_percent)
        one2many_indices = (D1byD2 > minimum_divide_percent) & (D2byD1.T > minimum_unique_percent)
        one2many_indices = one2many_indices & (one2many_indices.sum(axis = 1) >= 2).values[:, np.newaxis]
        many2one_indices = (D2byD1 > minimum_divide_percent) & (D1byD2.T > minimum_unique_percent)
        many2one_indices = many2one_indices & (many2one_indices.sum(axis = 1) >= 2).values[:, np.newaxis]
        row_index, col_index = np.where(indices & (indices.sum(axis = 1) == 1).values[:, np.newaxis] & (indices.sum(axis = 0) == 1).values & (~one2many_indices) & (~many2one_indices.T))
    #one2many
    if not D1byD2.empty:
        indices = (D1byD2 > minimum_divide_percent) & (D2byD1.T > minimum_unique_percent)
        row_index, col_index = np.where(indices & (indices.sum(axis = 1) >= 2).values[:, np.newaxis])
        while len(row_index) > 0:
            relation = pd.concat([relation, pd.DataFrame({'D1': D1byD2.index[row_index], 'relation': [ONE2MANY]*len(row_index), 'D2': D1byD2.columns[col_index]})])
            D1byD2, D2byD1 = _subset_and_normalize(D1byD2, D2byD1, relation, row_normalize)
            if D1byD2.empty:
                break
            indices = (D1byD2 > minimum_divide_percent) & (D2byD1.T > minimum_unique_percent)
            row_index, col_index = np.where(indices & (indices.sum(axis = 1) >= 2).values[:, np.newaxis])
    #many2one
    if not D1byD2.empty:
        indices = (D2byD1 > minimum_divide_percent) & (D1byD2.T > minimum_unique_percent)
        row_index, col_index = np.where(indices & (indices.sum(axis = 1) >= 2).values[:, np.newaxis])
        while len(row_index) > 0:
            relation = pd.concat([relation, pd.DataFrame({'D1': D2byD1.columns[col_index], 'relation': [MANY2ONE]*len(col_index), 'D2': D2byD1.index[row_index]})])
            D1byD2, D2byD1 = _subset_and_normalize(D1byD2, D2byD1, relation, row_normalize)
            if D1byD2.empty:
                break
            indices = (D2byD1 > minimum_divide_percent) & (D1byD2.T > minimum_unique_percent)
            row_index, col_index = np.where(indices & (indices.sum(axis = 1) >= 2).values[:, np.newaxis])
    #remaining
    if D1byD2.shape[0] > 0:
        relation = pd.concat([relation, pd.DataFrame({'D1': D1byD2.index, 'relation': [ONE2ONE]*len(D1byD2.index), 'D2': [REMAIN]*len(D1byD2.index)})])
    if D1byD2.shape[1] > 0:
        relation = pd.concat([relation, pd.DataFrame({'D1': [REMAIN]*len(D1byD2.columns), 'relation': [ONE2ONE]*len(D1byD2.columns), 'D2': D1byD2.columns})])
    return relation

def _reorder_dataset(sim_df: pd.DataFrame) -> np.ndarray:
    """
    For internal use. Reorder datasets based on their pairwise similarities.
    """
    stack_datasets = sim_df.iloc[sim_df.similarity.argmax(), [0, 1]].tolist()
    sim_df = pd.concat([sim_df, sim_df[['D2', 'D1', 'similarity']].rename(columns = {'D2':'D1', 'D1':'D2'})])
    remain_datasets = list(np.setdiff1d(np.unique(sim_df.D1), stack_datasets))
    while len(remain_datasets) >= 2:
        scores = [sim_df.loc[(sim_df.D1 == x) & sim_df.D2.isin(stack_datasets), 'similarity'].sum() for x in remain_datasets]
        stack_datasets.append(remain_datasets.pop(scores.index(max(scores))))
    stack_datasets.append(remain_datasets.pop())
    return np.array(stack_datasets, dtype = object)

class DistanceAlignment():
    """
    Class that performs cell type label harmonization across datasets.

    Parameters
    ----------
    base_distance
        A :class:`~celltypist.contro.distance.Distance` object.
    check
        Whether to check the supplied `base_distance` is correctly provided.
        (Default: `True`)
    dataset_order
        Order of datasets to be aligned. By default, the order is the same as that in the base distance matrix.
    row_normalize
        Whether to row normalize the confusion matrix to a sum of 1 in each iteration.
        (Default: `True`)
    minimum_unique_percent
        The minimum cell assignment fraction to claim a cell type as uniquely matched to a cell type from the other dataset.
        (Default: `0.5`)
    minimum_divide_percent
        The minimum cell assignment fraction to claim a cell type as divisible into two or more cell types from the other dataset.
        (Default: `0.1`)
    maximum_novel_percent
        The maximum cell assignment fraction to claim a cell type as novel to a given dataset.
        (Default: `0.05`)

    Attributes
    ----------
    base_distance
        The :class:`~celltypist.contro.distance.Distance` object.
    dataset_order
        Order of datasets to be aligned.
    row_normalize
        Whether to row normalize the confusion matrix to a sum of 1 in each iteration.
    minimum_unique_percent
        The minimum cell assignment fraction to claim a cell type as uniquely matched to a cell type from the other dataset.
    minimum_divide_percent
        The minimum cell assignment fraction to claim a cell type as divisible into two or more cell types from the other dataset.
    maximum_novel_percent
        The maximum cell assignment fraction to claim a cell type as novel to a given dataset.
    relation
        A :class:`~pandas.DataFrame` representing the harmonization result.
    aligned_datasets
        List of datasets that are already harmonized.
    groups
        Cell type groups (high-hierarchy cell types) categorizing the rows of `.relation`.
    reannotation
        A :class:`~pandas.DataFrame` representing the reannotated cell types.
    minimum_unique_percents
        List of `minimum_unique_percent` values which are used along harmonization iterations in order to get the best alignment.
        This attribute is obtained through the :func:`~celltypist.contro.align.DistanceAlignment.best_align` function.
    minimum_divide_percents
        List of `minimum_divide_percent` values which are used along harmonization iterations in order to get the best alignment.
        This attribute is obtained through the :func:`~celltypist.contro.align.DistanceAlignment.best_align` function.
    """
    def __init__(self, base_distance: Distance, check: bool = True, dataset_order: Optional[Union[list, tuple, np.ndarray, pd.Series, pd.Index]] = None,
                 row_normalize: bool = True, minimum_unique_percent: float = 0.5, minimum_divide_percent: float = 0.1, maximum_novel_percent: float = 0.05):
        if check:
            if not isinstance(base_distance, Distance) or not base_distance.symmetric():
                raise ValueError(
                        f"🛑 Please provide a symmetric `Distance` object")
            if not hasattr(base_distance, 'assignment'):
                raise AttributeError(
                        f"🛑 No `.assignment` attribute in the `base_distance`. Apply the `.assign` method first")
        self.base_distance = base_distance
        if dataset_order is None:
            dataset_order = np.unique(self.base_distance.cell_type.dataset)
        else:
            dataset_order = np.array(dataset_order, dtype = object)
            if not np.array_equal(np.sort(dataset_order), np.unique(self.base_distance.cell_type.dataset)):
                raise ValueError(
                        f"🛑 Please provide a comprehensive order of datasets with correct names")
        self.dataset_order = dataset_order
        self.row_normalize = row_normalize
        self.minimum_unique_percent = minimum_unique_percent
        self.minimum_divide_percent = minimum_divide_percent
        self.maximum_novel_percent = maximum_novel_percent

    def reorder_dataset(self, weights: Union[list, tuple, np.ndarray, pd.Series, pd.Index] = (2, 1, -1, -2), return_similarity: bool = False) -> Union[None, pd.DataFrame]:
        """
        Reorder the datasets such that similar datasets will be harmonized first. This method can also be used to calculate CellTypist-defined inter-dataset similarities.

        Parameters
        ----------
        weights
            Weights assigned to one-to-one, one/many-to-many/one, novel, and remaining cell type matches, respectively. Default to 2, 1, -1, -2.
            Inter-cell-type similarities will be weighted by these values to derive the weighted sum of similarity between each pair of datasets.
        return_similarity
            Whether to return the data frame of dataset-dataset similarities.
            (Default: `False`)

        Returns
        ----------
        Reordered datasets as the attribute `.dataset_order` and if `return_similarity = True`, return a :class:`~pandas.DataFrame` of dataset-dataset similarities.
        """
        if len(self.dataset_order) == 2:
            logger.warn(f"⚠️ Warning: only two datasets exist, no need to reorder them")
            return
        weights = np.array(weights, dtype = 'float')
        meta_distance = self.base_distance.to_meta(False)
        meta_similarity = 1 - meta_distance
        sim_df = pd.DataFrame([self.dataset_order[[i, j]] for i in range(0, len(self.dataset_order) - 1) for j in range(i + 1, len(self.dataset_order))], columns = ['D1', 'D2'])
        scores = []
        for _, s in sim_df.iterrows():
            D1 = s.values[0]
            D2 = s.values[1]
            sub_meta_similarity = meta_similarity.loc[meta_similarity.index.str.startswith(D1 + ": "), meta_similarity.columns.str.startswith(D2 + ": ")]
            sub_meta_similarity.index = sub_meta_similarity.index.str.replace(D1 + ": ", '', regex = False)
            sub_meta_similarity.columns = sub_meta_similarity.columns.str.replace(D2 + ": ", '', regex = False)
            relation = self.pairwise_align(D1, D2, False)
            ss = np.full(relation.shape[0], -1, dtype = 'float')
            ws = np.full(relation.shape[0], weights[0], dtype = 'float')
            #REMAIN NOVEL
            flag_REMAIN_1 = (relation[D1] == REMAIN).values
            flag_REMAIN_2 = (relation[D2] == REMAIN).values
            flag_NOVEL_1 = (relation[D1] == NOVEL).values
            flag_NOVEL_2 = (relation[D2] == NOVEL).values
            flag_NONE_1 = flag_REMAIN_1 | flag_NOVEL_1
            if flag_NONE_1.sum() > 0:
                ss[flag_NONE_1] = 1 - sub_meta_similarity[relation.loc[flag_NONE_1, D2].values].max(axis = 0).values
            flag_NONE_2 = flag_REMAIN_2 | flag_NOVEL_2
            if flag_NONE_2.sum() > 0:
                ss[flag_NONE_2] = 1 - sub_meta_similarity.loc[relation.loc[flag_NONE_2, D1].values].max(axis = 1).values
            ws[flag_REMAIN_1 | flag_REMAIN_2] = weights[3]
            ws[flag_NOVEL_1 | flag_NOVEL_2] = weights[2]
            #other ss
            other_flag = ~(flag_NONE_1 | flag_NONE_2)
            if other_flag.sum() > 0:
                ss[other_flag] = np.diag(sub_meta_similarity.loc[relation.loc[other_flag, D1].values, relation.loc[other_flag, D2].values])
            ##remove in the future-->
            assert np.all(ss != -1)
            ##<<-remove in the future
            #ONE2MANY and MANY2ONE weights
            flag_ONE2MANY = (relation.relation == ONE2MANY).values
            if flag_ONE2MANY.sum() > 0:
                DS = relation.loc[flag_ONE2MANY, D1].value_counts()
                ws[flag_ONE2MANY] = (weights[1] / DS)[relation.loc[flag_ONE2MANY, D1].values].values
            flag_MANY2ONE = (relation.relation == MANY2ONE).values
            if flag_MANY2ONE.sum() > 0:
                DS = relation.loc[flag_MANY2ONE, D2].value_counts()
                ws[flag_MANY2ONE] = (weights[1] / DS)[relation.loc[flag_MANY2ONE, D2].values].values
            #final scores
            scores.append(np.sum(ss * ws) / sum(sub_meta_similarity.shape))
        sim_df['similarity'] = scores
        self.dataset_order = _reorder_dataset(sim_df)
        if return_similarity:
            return sim_df

    def pairwise_align(self, D1: str, D2: str, check: bool = True) -> pd.DataFrame:
        """
        Pairwise alignment of cell types between two datasets.

        Parameters
        ----------
        D1
            Name of the first dataset.
        D2
            Name of the second dataset.
        check
            Whether to check names of the two datasets are contained in the :attr:`~celltypist.contro.align.DistanceAlignment.base_distance`.
            (Default: `True`)

        Returns
        ----------
        :class:`~pandas.DataFrame`
            A :class:`~pandas.DataFrame` with three columns:
            1) **name of dataset 1**, cell types from dataset 1.
            2) **relation**, being either '=', '∋' or '∈'.
            3) **name of dataset 2**, cell types from dataset 2.
        """
        D1byD2, D2byD1 = self.base_distance.to_pairwise_confusion(D1, D2, check)
        relation = _pairwise_align(D1byD2, D2byD1, False, self.row_normalize, self.minimum_unique_percent, self.minimum_divide_percent, self.maximum_novel_percent)
        return relation.rename(columns = {'D1': D1, 'D2': D2})

    def multi_align(self, relation: pd.DataFrame, D: str, check: bool = True) -> pd.DataFrame:
        """
        Multiple alignment of cell types across datasets. Cell types from a new dataset will be integrated into the previous harmonization data frame.

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
        :class:`~pandas.DataFrame`
            A :class:`~pandas.DataFrame` with multiple columns:
            1) **name of dataset 1**, cell types from dataset 1.
            2) **relation**, being either '=', '∋' or '∈'.
            3) **name of dataset 2**, cell types from dataset 2.
            4) ...
            N) **name of the new dataset**, cell types from the new dataset.
        """
        D1byD2, D2byD1 = self.base_distance.to_multi_confusion(relation, D, check)
        new_relation = _pairwise_align(D1byD2, D2byD1, False, self.row_normalize, self.minimum_unique_percent, self.minimum_divide_percent, self.maximum_novel_percent)
        relation.index = relation.apply(lambda row: ' '.join(row.values), axis = 1).values
        relation.loc[NOVEL] = np.tile([NOVEL, ONE2ONE], int(relation.shape[1]/2)+1)[:-1]
        relation.loc[REMAIN] = np.tile([REMAIN, ONE2ONE], int(relation.shape[1]/2)+1)[:-1]
        relation = relation.loc[new_relation.D1.values]
        relation.index = new_relation.index
        relation['__relation'] = new_relation['relation'].values
        relation[D] = new_relation['D2'].values
        return relation.rename(columns = {'__relation': 'relation'})

    def align(self, datasets: Optional[Union[list, tuple, np.ndarray, pd.Series, pd.Index]] = None) -> None:
        """
        Iterative alignment of cell types across datasets.

        Parameters
        ----------
        datasets
            Datasets to be aligned. Default to using all datasets available.

        Returns
        ----------
        None
            A :class:`~pandas.DataFrame` with multiple columns added as the attribute `.relation`:
            1) **name of dataset 1**, cell types from dataset 1.
            2) **relation**, being either '=', '∋' or '∈'.
            3) **name of dataset 2**, cell types from dataset 2.
            4) ...
            N) **name of the last dataset**, cell types from the last dataset.
        """
        if datasets is None:
            datasets = self.dataset_order
        else:
            datasets = np.array(datasets, dtype = object)
            if not set(datasets).issubset(self.dataset_order):
                raise ValueError(
                        f"🛑 The following datasets are not found: {set(datasets).difference(self.dataset_order)}")
        logger.info(f"🖇 Harmonizing cell types of {datasets[0]} and {datasets[1]}")
        relation = self.pairwise_align(datasets[0], datasets[1], False)
        if len(datasets) >= 3:
            for i in range(2, len(datasets)):
                logger.info(f"🖇 Harmonizing cell types of {datasets[i]}")
                relation = self.multi_align(relation, datasets[i], False)
        logger.info(f"✅ Harmonization done!")
        self.relation = relation

    @property
    def aligned_datasets(self) -> np.ndarray:
        """Get the datasets which are already harmonized."""
        return np.array(self.relation.columns[0::2])

    @property
    def groups(self) -> np.ndarray:
        """Get the cell type groups (high hierarchy) based on the relation table."""
        return _identify_relation_groups(self.relation, order_row = False, order_column = False)[0]

    def __repr__(self):
        base = f"Cross-dataset cell type alignment for {len(self.dataset_order)} datasets"
        base += f"\n    base_distance: a cross-dataset distance object"
        if hasattr(self, 'relation'):
            base += f"\n    aligned_datasets: {str(list(self.aligned_datasets))[1:-1]}"
            base += f"\n    relation: data frame of the harmonization table"
        if hasattr(self, 'reannotation'):
            base += f"\n    reannotation: data frame of the reannotated cells ({str(list(self.reannotation.columns))[1:-1]})"
        return base

    def update(self, datasets: Optional[Union[str, list, tuple, np.ndarray, pd.Series, pd.Index]] = None) -> None:
        """
        Iteratively update the alignment of cell types across datasets.

        Parameters
        ----------
        datasets
            Datasets to be aligned. Default to using all the remaining datasets.

        Returns
        ----------
        None
            An updated :class:`~pandas.DataFrame` with multiple columns added as the attribute `.relation`:
            1) **name of dataset 1**, cell types from dataset 1.
            2) **relation**, being either '=', '∋' or '∈'.
            3) **name of dataset 2**, cell types from dataset 2.
            4) ...
            N) **name of the last dataset**, cell types from the last dataset.
        """
        if not hasattr(self, 'relation'):
            raise AttributeError(
                    f"🛑 No harmonization result exists. Please run the `.align` method first")
        if len(self.aligned_datasets) == len(self.dataset_order):
            logger.warn(f"⚠️ All datasets have been harmonized. No update is needed")
            return
        remaining_datasets = self.dataset_order[~np.isin(self.dataset_order, self.aligned_datasets)]
        if datasets is None:
            datasets = remaining_datasets
        else:
            datasets = [datasets] if isinstance(datasets, str) else np.array(datasets, dtype = object)
            if not set(datasets).issubset(remaining_datasets):
                raise ValueError(
                        f"🛑 Please provide dataset names from the following list: {set(remaining_datasets)}")
        for dataset in datasets:
            logger.info(f"🖇 Harmonizing cell types of {dataset}")
            self.relation = self.multi_align(self.relation, dataset, False)
        logger.info(f"✅ Harmonization done!")

    def best_align(self, dataset_order: Optional[Union[list, tuple, np.ndarray, pd.Series, pd.Index]] = None, minimum_unique_percents: Union[list, tuple, np.ndarray, pd.Series, pd.Index, float] = (0.4, 0.5, 0.6, 0.7, 0.8), minimum_divide_percents: Union[list, tuple, np.ndarray, pd.Series, pd.Index, float] = (0.1, 0.15, 0.2)):
        """
        Iterative alignment of cell types across datasets by finding the best parameter combo in each iteration.

        Parameters
        ----------
        dataset_order
            Order of datasets to be aligned. This can also be a subset of datasets.
            Default to the dataset order in the `DistanceAlignment` object.
        minimum_unique_percents
            The minimum cell assignment fraction(s) to claim a cell type as uniquely matched to a cell type from the other dataset.
            By default, five values will be tried (0.4, 0.5, 0.6, 0.7, 0.8) to find the one that produces least alignments in each harmonization iteration.
        minimum_divide_percents
            The minimum cell assignment fraction(s) to claim a cell type as divisible into two or more cell types from the other dataset.
            By default, three values will be tried (0.1, 0.15, 0.2) to find the one that produces least alignments in each harmonization iteration.

        Returns
        ----------
        None
            A :class:`~pandas.DataFrame` with multiple columns added as the attribute `.relation`:
            1) **name of dataset 1**, cell types from dataset 1.
            2) **relation**, being either '=', '∋' or '∈'.
            3) **name of dataset 2**, cell types from dataset 2.
            4) ...
            N) **name of the last dataset**, cell types from the last dataset.
        """
        if dataset_order is not None:
            dataset_order = np.array(dataset_order, dtype = object)
            if not set(dataset_order).issubset(self.dataset_order):
                raise ValueError(
                        f"🛑 The following dataset(s) are not found: {set(dataset_order).difference(self.dataset_order)}")
        else:
            dataset_order = self.dataset_order
        original_mup = self.minimum_unique_percent
        original_mdp = self.minimum_divide_percent
        minimum_unique_percents = np.array([minimum_unique_percents]) if isinstance(minimum_unique_percents, float) else np.array(minimum_unique_percents)
        minimum_divide_percents = np.array([minimum_divide_percents]) if isinstance(minimum_divide_percents, float) else np.array(minimum_divide_percents)
        mups = np.full(len(dataset_order) - 1, -1, dtype = 'float')
        mdps = mups.copy()
        logger.info(f"🖇 Harmonizing cell types of {dataset_order[0]} and {dataset_order[1]}")
        n_rows = np.inf
        for minimum_unique_percent in minimum_unique_percents:
            for minimum_divide_percent in minimum_divide_percents:
                self.minimum_unique_percent = minimum_unique_percent
                self.minimum_divide_percent = minimum_divide_percent
                per_relation = self.pairwise_align(dataset_order[0], dataset_order[1], check = False)
                if per_relation.shape[0] <= n_rows:
                    relation = per_relation
                    n_rows = per_relation.shape[0]
                    mups[0] = minimum_unique_percent
                    mdps[0] = minimum_divide_percent
        if len(dataset_order) >= 3:
            for i in range(2, len(dataset_order)):
                logger.info(f"🖇 Harmonizing cell types of {dataset_order[i]}")
                n_rows = np.inf
                for minimum_unique_percent in minimum_unique_percents:
                    for minimum_divide_percent in minimum_divide_percents:
                        self.minimum_unique_percent = minimum_unique_percent
                        self.minimum_divide_percent = minimum_divide_percent
                        per_relation = self.multi_align(relation.copy(), dataset_order[i], check = False)
                        if per_relation.shape[0] <= n_rows:
                            expand_relation = per_relation
                            n_rows = per_relation.shape[0]
                            mups[i - 1] = minimum_unique_percent
                            mdps[i - 1] = minimum_divide_percent
                relation = expand_relation
        self.minimum_unique_percents = mups
        self.minimum_divide_percents = mdps
        self.minimum_unique_percent = original_mup
        self.minimum_divide_percent = original_mdp
        self.relation = relation

    def reannotate(self, show_iteration: bool = False, add_group: bool = True, prefix: str = '') -> None:
        """
        Reannotate each cell into the harmonized cell type.

        Parameters
        ----------
        show_iteration
            Whether to store the cell type reannotation result for each harmonization iteration.
            (Default: `False`)
        add_group
            Whether to annotate out cell type group information.
            (Default: `True`)
        prefix
            Prefix of the harmonization columns for all iterations. Default to no prefix.

        Returns
        ----------
        None
            A :class:`~pandas.DataFrame` with multiple columns added as the attribute `.reannotation`:
            1) **dataset**, datasets where the cells are from.
            2) **cell_type**, cell types annotated by the original studies/datasets.
            3) **roundN** or **reannotation**, prefixed with `prefix`; cell types reannotated by the harmonization process.
            4) **group**, prefixed with `prefix`; annotated cell type groups.
        """
        if not hasattr(self, 'relation'):
            raise AttributeError(
                    f"🛑 No harmonization result (`.relation`) exists")
        reannotation = self.base_distance.cell.set_index('ID', inplace = False, drop = True)
        lend = len(self.aligned_datasets)
        reannotation[[f"{prefix}round{x}" for x in range(1, lend)]] = UNASSIGN
        assignment = self.base_distance.assignment[self.aligned_datasets]
        for _, s in self.relation.iterrows():
            celltypes = s.values[0::2]
            flags = (assignment == celltypes) | np.isin(celltypes, [NOVEL, REMAIN])
            non_existing_datasets = self.aligned_datasets[np.isin(celltypes, [NOVEL, REMAIN])]
            flags.loc[reannotation.dataset.isin(non_existing_datasets).values, :] = False
            for N in range(1, lend):
                if (not show_iteration) and (N != lend - 1):
                    continue
                sub_s = s.values[:2*N+1]
                sub_celltypes = sub_s[0::2]
                if np.all(sub_celltypes == NOVEL) or np.all(sub_celltypes == REMAIN):
                    continue
                reannotation.loc[np.all(flags[self.aligned_datasets[:N+1]], axis = 1).values, f"{prefix}round{N}"] = ' '.join(sub_s)
        if add_group:
            groups, new_relation = _identify_relation_groups(self.relation, order_row = False, order_column = False)
            group_mapping = dict()
            for i in range(len(groups)):
                group_mapping.update({j: groups[i] for j in new_relation.iloc[i].values})
            reannotation[f"{prefix}group"] = (reannotation.dataset + SEP1 + reannotation.cell_type).replace(group_mapping)
            reannotation.loc[reannotation[f"{prefix}group"].str.contains(SEP1), f"{prefix}group"] = UNASSIGN
        if not show_iteration:
            reannotation.rename(columns = {f"{prefix}round{lend-1}": f"{prefix}reannotation"}, inplace = True)
            if lend >= 3:
                reannotation.drop(columns = [f"{prefix}round{x}" for x in range(1, lend-1)], inplace = True)
        self.reannotation = reannotation

    @staticmethod
    def load(alignment_file: str):
        """Load the DistanceAlignment file."""
        with open(alignment_file, "rb") as fh:
            return pickle.load(fh)

    def write(self, file: str) -> None:
        """Write out the DistanceAlignment."""
        file = os.path.splitext(file)[0] + '.pkl'
        with open(file, 'wb') as output:
            pickle.dump(self, output)
