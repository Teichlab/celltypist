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

def _expand_ancestor_sets(y_true: np.ndarray, y_pred: np.ndarray, tree: Tree, include_root: bool) -> tuple:
    """
    For internal use. Expand true and predicted labels into ancestor sets for each sample.
    """
    label_path = {}
    for label in np.unique(np.concatenate([y_true, y_pred])):
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
                f"🛑 `metric_type` must be one of `'accuracy'`, `'precision'`, `'recall'`, or `'f1'`")
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

def _compute_micro(true_sets: list, pred_sets: list, metric_type: str) -> float:
    """
    For internal use. Compute hierarchical micro-metrics via global aggregation over all samples.
    """
    if metric_type not in ('accuracy', 'precision', 'recall', 'f1'):
        raise ValueError(
                f"🛑 `metric_type` must be one of `'accuracy'`, `'precision'`, `'recall'`, or `'f1'`")
    total_intersection = 0
    total_true = 0
    total_pred = 0
    total_union = 0
    for ts, ps in zip(true_sets, pred_sets):
        total_intersection += len(ts & ps)
        total_true += len(ts)
        total_pred += len(ps)
        total_union += len(ts | ps)
    if metric_type == "accuracy":
        return total_intersection / total_union if total_union else 0.0
    elif metric_type == "precision":
        return total_intersection / total_pred if total_pred else 0.0
    elif metric_type == "recall":
        return total_intersection / total_true if total_true else 0.0
    else:
        denominator = (total_true + total_pred) / 2
        return total_intersection / denominator if denominator else 0.0

def hier_accuracy(y_true: Union[list, tuple, np.ndarray, pd.Series, pd.Index], y_pred: Union[list, tuple, np.ndarray, pd.Series, pd.Index], tree: Tree, include_root: bool = False,
                  average: Optional[str] = 'sample macro') -> Union[tuple, float]:
    """
    Compute hierarchical accuracy based on Jaccard similarity between the expanded ancestor sets of true and predicted labels.

    Parameters
    ----------
    y_true
        Ground-truth labels.
    y_pred
        Predicted labels.
    tree
        A :class:`~celltypist.tree.Tree` object representing the predefined cell type hierarchy.
    include_root
        Whether to include the root node when constructing ancestor sets.
        (Default: `False`)
    average
        Averaging strategy:
        1) 'micro': compute the global Jaccard score by aggregating intersections and unions across all samples.
        2) 'sample macro': compute the mean of per-sample Jaccard scores (default).
        3) 'label macro': compute the mean of per-label averaged scores.
        4) None: return per-label hierarchical accuracy values.

    Returns
    ----------
    Union[tuple, float]
        Returns a tuple containing labels and their corresponding hierarchical accuracy values (`average = None`). Otherwise, returns a single aggregated hierarchical accuracy.
    """
    y_true, y_pred = _check_tree_and_labels(y_true, y_pred, tree)
    if average not in (None, 'sample macro', 'label macro', 'micro'):
        raise ValueError(
                f"🛑 `average` must be one of `None`, `'sample macro'`, `'label macro'`, or `'micro'`")
    true_sets, pred_sets = _expand_ancestor_sets(y_true, y_pred, tree, include_root)
    if average == 'micro':
        return _compute_micro(true_sets, pred_sets, 'accuracy')
    else:
        return _compute_macro(true_sets, pred_sets, 'accuracy', average, y_true)
