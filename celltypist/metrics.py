import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score
from typing import Union
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
