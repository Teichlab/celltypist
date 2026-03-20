import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
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
    if labels is None:
        labels = np.unique(y_true)
    else:
        labels = np.asarray(labels)
        diff = set(labels).difference(tree.cell_types(leaf_only = False))
        if diff:
            raise ValueError(
                    f"🛑 The following labels are not in the tree: {diff}")
    scores = precision_score(y_true, y_pred, labels = labels, average = average, sample_weight = None, zero_division = 0)
    if average is None:
        return labels, scores
    else:
        return scores

def flat_recall(y_true: Union[list, tuple, np.ndarray, pd.Series, pd.Index], y_pred: Union[list, tuple, np.ndarray, pd.Series, pd.Index], tree: Tree,
                labels: Optional[Union[list, tuple, np.ndarray, pd.Series, pd.Index]] = None, average: Optional[str] = None) -> Union[tuple, float]:
    """
    Compute standard flat classification recall, ignoring hierarchical relationships between cell type labels.

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
        Strategies for recall calculation: `micro` computes recall globally by aggregating TP and FN from all labels, whereas `macro` computes the mean of per-label recall values.
        If no averaging strategy is specified, recall is calculated individually for each label in `labels`.

    Returns
    ----------
    Union[tuple, float]
        Returns either a tuple containing labels and their corresponding recall values (default), or a single recall value when `average = 'micro'` or `average = 'macro'`.
    """
    y_true, y_pred = _check_tree_and_labels(y_true, y_pred, tree)
    if average not in (None, 'micro', 'macro'):
        raise ValueError(
                f"🛑 If specified, `average` must be either `'micro'` or `'macro'`")
    if labels is None:
        labels = np.unique(y_true)
    else:
        labels = np.asarray(labels)
        diff = set(labels).difference(tree.cell_types(leaf_only = False))
        if diff:
            raise ValueError(
                    f"🛑 The following labels are not in the tree: {diff}")
    scores = recall_score(y_true, y_pred, labels = labels, average = average, sample_weight = None, zero_division = 0)
    if average is None:
        return labels, scores
    else:
        return scores

def flat_f1(y_true: Union[list, tuple, np.ndarray, pd.Series, pd.Index], y_pred: Union[list, tuple, np.ndarray, pd.Series, pd.Index], tree: Tree,
            labels: Optional[Union[list, tuple, np.ndarray, pd.Series, pd.Index]] = None, average: Optional[str] = None) -> Union[tuple, float]:
    """
    Compute standard flat classification F1 score, ignoring hierarchical relationships between cell type labels.

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
        Strategies for F1 score calculation: `micro` computes F1 score globally by aggregating TP, FP, and FN from all labels, whereas `macro` computes the mean of per-label F1 scores.
        If no averaging strategy is specified, F1 score is calculated individually for each label in `labels`.

    Returns
    ----------
    Union[tuple, float]
        Returns either a tuple containing labels and their corresponding F1 scores (default), or a single F1 score when `average = 'micro'` or `average = 'macro'`.
    """
    y_true, y_pred = _check_tree_and_labels(y_true, y_pred, tree)
    if average not in (None, 'micro', 'macro'):
        raise ValueError(
                f"🛑 If specified, `average` must be either `'micro'` or `'macro'`")
    if labels is None:
        labels = np.unique(y_true)
    else:
        labels = np.asarray(labels)
        diff = set(labels).difference(tree.cell_types(leaf_only = False))
        if diff:
            raise ValueError(
                    f"🛑 The following labels are not in the tree: {diff}")
    scores = f1_score(y_true, y_pred, labels = labels, average = average, sample_weight = None, zero_division = 0)
    if average is None:
        return labels, scores
    else:
        return scores

def _expand_ancestor_sets(y_true: np.ndarray, y_pred: np.ndarray, tree: Tree, labels: np.ndarray, include_root: bool) -> tuple:
    """
    For internal use. Expand true and predicted labels into ancestor sets for each sample.
    """
    label_path = {}
    for label in labels:
        node_path = tree.extract_path(label, print_path = False)
        names = {node.original_name for node in node_path}
        if not include_root:
            names.discard(tree.root.original_name)
        label_path[label] = names
    return [label_path[t] for t in y_true], [label_path[p] for p in y_pred]

def _compute_macro(true_sets: list, pred_sets: list, metric_type: str, average: Union[str, None], y_true: Optional[np.ndarray] = None) -> Union[tuple, float]:
    """
    For internal use. Compute hierarchical macro-metrics with different averaging strategies.
    """
    if metric_type not in ('accuracy', 'precision', 'recall', 'f1'):
        raise ValueError(
                f"🛑 `metric_type` must be one of `'accuracy'`, `'precision'`, `'recall'`, and `'f1'`")
    if average not in (None, 'sample macro', 'label macro'):
        raise ValueError(
                f"🛑 If specified, `average` must be either `'sample macro'` or `'label macro'`")
    scores = []
    for ts, ps in zip(true_sets, pred_sets):
        if metric_type == 'accuracy':
            denominator = len(ts | ps)
        elif metric_type == 'precision':
            denominator = len(ps)
        elif metric_type == 'recall':
            denominator = len(ts)
        else:
            denominator = (len(ts) + len(ps)) / 2
        scores.append(len(ts & ps) / denominator if denominator else 0.0)
    scores = np.asarray(scores)
    if average == 'sample macro':
        return float(np.mean(scores))
    else:
        labels = np.unique(y_true)
        values = np.asarray([float(np.mean(scores[y_true == label])) for label in labels])
        if average is None:
            return labels, values
        else:
            return float(np.mean(values))

#def hier_accuracy(y_true: Union[list, tuple, np.ndarray, pd.Series, pd.Index], y_pred: Union[list, tuple, np.ndarray, pd.Series, pd.Index], tree: Tree, include_root: bool = False,
#                  average: Optional[str] = 'sample macro') -> Union[tuple, float]:
#    """
#    Compute hierarchical accuracy based on Jaccard similarity between expanded ancestor sets of true and predicted labels.
#
#    Parameters
#    ----------
#    y_true
#        Ground-truth labels.
#    y_pred
#        Predicted labels.
#    tree
#        A :class:`~celltypist.tree.Tree` object representing the predefined cell type hierarchy.
#    include_root
#        Whether to include the root node when constructing ancestor sets.
#        (Default: `False`)
#    average
#        Averaging strategy:
#        - 'micro': 
#        - 'sample macro': mean of per-sample Jaccard scores (default)
#        - 'label macro': mean of per-label averaged scores
#        - None: return per-label scores
#        (Default: `'sample macro'`)
#
#    Returns
#    ----------
#    Union[tuple, float]
#        Returns a tuple containing labels and their corresponding hierarchical accuracy values (`average = None`). Otherwise returns a single aggregated hierarchical accuracy.
#    """
#    y_true, y_pred = _check_tree_and_labels(y_true, y_pred, tree)
#    true_sets, pred_sets = _expand_ancestor_sets(y_true, y_pred, tree, np.unique(np.concatenate([y_true, y_pred])), include_root)
#    return _compute_macro(true_sets, pred_sets, 'accuracy', average, y_true)
