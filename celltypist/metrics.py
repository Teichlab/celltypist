import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, precision_score
from typing import Union, Optional
from .tree import Tree

def _check_tree_and_labels(y_true: Union[list, tuple, np.ndarray, pd.Series, pd.Index], y_pred: Union[list, tuple, np.ndarray, pd.Series, pd.Index], tree: Tree) -> tuple:
    """
    For internal use. Validate that inputs are compatible with the provided tree.
    """
    if not isinstance(tree, Tree):
        raise TypeError(
                f"🛑 Please provide a `Tree` instance for the `tree` argument")
    all_labels = tree.cell_types(leaf_only = False)
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    diff_true = set(y_true).difference(all_labels)
    if diff_true:
        raise ValueError(
                f"🛑 The following ground-truth labels are not in the tree: {diff_true}")
    diff_pred = set(y_pred).difference(all_labels)
    if diff_pred:
        raise ValueError(
                f"🛑 The following predicted labels are not in the tree: {diff_pred}")
    return y_true, y_pred

def flat_accuracy(y_true: Union[list, tuple, np.ndarray, pd.Series, pd.Index], y_pred: Union[list, tuple, np.ndarray, pd.Series, pd.Index], tree: Tree) -> float:
    """
    Compute standard flat classification accuracy, ignoring hierarchical relationships between cell type labels.

    Parameters
    ----------
    y_true
        Ground-truth labels.
    y_pred
        Predicted labels.
    tree
        A :class:`~celltypist.tree.Tree` object representing the predefined cell type hierarchy.

    Returns
    ----------
    float
        Fraction of correctly classified labels.
    """
    y_true, y_pred = _check_tree_and_labels(y_true, y_pred, tree)
    return accuracy_score(y_true, y_pred, normalize = True, sample_weight = None)

def flat_precision(y_true: Union[list, tuple, np.ndarray, pd.Series, pd.Index], y_pred: Union[list, tuple, np.ndarray, pd.Series, pd.Index], tree: Tree,
                   labels: Optional[Union[list, tuple, np.ndarray, pd.Series, pd.Index]] = None, average: Optional[str] = None) -> Union[tuple, float]:
    """
    Compute standard flat classification precision, ignoring hierarchical relationships between cell type labels.

    Parameters
    ----------
    y_true
        Ground-truth labels.
    y_pred
        Predicted labels.
    tree
        A :class:`~celltypist.tree.Tree` object representing the predefined cell type hierarchy.
    labels
        The set of labels to include in the evaluation. Defaults to the unique labels present in `y_true`.
    average
        Strategies for precision calculation: `micro` computes precision globally by aggregating TP and FP from all labels, whereas `macro` computes the mean of per-label precision values.
        If no averaging strategy is specified, precision is calculated individually for each label in `labels`.

    Returns
    ----------
    Union[tuple, float]
        Returns either a tuple containing labels and their corresponding precision values (default), or a single precision value when `average = 'micro'` or `average = 'macro'`.
    """
    y_true, y_pred = _check_tree_and_labels(y_true, y_pred, tree)
    if average not in (None, 'micro', 'macro'):
        raise ValueError(
                f"🛑 If specified, `average` must be either `'micro'` or `'macro'`")
    labels = np.asarray(labels) if labels is not None else np.unique(y_true)
    scores = precision_score(y_true, y_pred, labels = labels, average = average, sample_weight = None, zero_division = 0)
    if average is None:
        return labels, scores
    else:
        return scores
