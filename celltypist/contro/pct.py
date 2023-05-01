from sklearn.tree._tree import TREE_LEAF, TREE_UNDEFINED
from sklearn.tree import DecisionTreeRegressor
from sklearn import __version__ as skv
import numpy as np
from scipy.stats import f
from typing import Optional, Union
from .. import logger

class PredictiveClusteringTree(DecisionTreeRegressor):
    """
    Class that uses predictive clustering trees (PCT) for multi-output prediction.
    Note this is a specialized PCT with the prototype as the mean vector and the distance measure as the (sum of) intra-cluster variance.
    Such a PCT is equivalent to CART regression tree (with specialized parameters and post-pruning).

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

    Attributes
    ----------
    n_features_in_
        Number of features.
    n_outputs_
        Number of outputs.
    tree_
        A :class:`~sklearn.tree._tree.Tree` object structured by parallel arrays.
    p_value
        F-test-based p-value for each node or leaf in the tree.
    """
    def __init__(self,
                 *,
                 max_depth: Optional[int] = None,
                 min_samples_split: Union[int, float] = 20,
                 min_samples_leaf: Union[int, float] = 10,
                 min_weight_fraction_leaf: float = 0.0,
                 random_state: Optional[int] = None,
                 max_leaf_nodes: Optional[int] = None,
                 F_test_prune: bool = True,
                 p_thres: float = 0.05):
        criterion = 'mse' if skv.split('.')[0] == '0' else 'squared_error'
        super().__init__(criterion = criterion,                             #L2 norm multi-output decision tree = pct with mean as prototype and SSe as distance measure
                         splitter = 'best',
                         max_depth = max_depth,
                         min_samples_split = min_samples_split,             #min_samples_leaf*2
                         min_samples_leaf = min_samples_leaf,               #minimum cell cluster size
                         min_weight_fraction_leaf = min_weight_fraction_leaf,
                         max_features = None,                               #use all features
                         random_state = random_state,
                         max_leaf_nodes = max_leaf_nodes,
                         min_impurity_decrease = 0.0,                       #split as much as possible -> no good way to select a proper value
                         ccp_alpha = 0.0)                                   #use F test for pruning, ignore ccp pruning
        self.F_test_prune = F_test_prune
        self.p_thres = p_thres

    def fit(self, X, y, sample_weight = None) -> None:
        """
        Fit a PCT with the training dataset.

        Parameters
        ----------
        X
            Sample-by-feature array-like matrix.
        y
            Sample-by-output array-like matrix.
        sample_weight
            Sample weights. Default to equal sample weights.

        Returns
        -------
        None
            Fitted and (possibly) pruned tree.
        """
        if isinstance(X, spmatrix) and ((X.indices.dtype == 'int64') or (X.indptr.dtype == 'int64')):
            super().fit(X.toarray(), y, sample_weight = sample_weight, check_input = True)
        else:
            super().fit(X, y, sample_weight = sample_weight, check_input = True)
        if self.F_test_prune:
            self.F_test()
            self.prune_tree(p_thres = self.p_thres)

    def is_node(self, index: int) -> bool:
        """
        Check whether a given index is a node.

        Parameters
        ----------
        index
            Index of the node/leaf in the arrays of the tree structure.

        Returns
        -------
        bool
            `True` or `False` indicating whether the given index is an internal node.
        """
        return self.tree_.children_left[index] != self.tree_.children_right[index]

    def F_test(self) -> None:
        """
        F test for each internal node.
        For each node, the corresponding F distribution has the degrees of freedom `n_output * (n_sample - 1)` and `n_output * (n_sample - 2)`, and the value (q) of `node_impurity * n_sample / (n_sample - 1)` divided by `(left_child_impurity * left_n_sample + right_child_impurity * right_n_sample) / (n_sample - 2)`.

        Returns
        -------
        None
            Modified tree with F-test p-values. Leaves are assigned 1 constantly.
        """
        ps = np.ones(self.tree_.node_count)
        for x in range(self.tree_.node_count):
            if self.is_node(x):
                f_score_numerator = self.tree_.impurity[x] * self.tree_.n_node_samples[x] / (self.tree_.n_node_samples[x] - 1)
                f_score_denominator = (self.tree_.impurity[self.tree_.children_left[x]]*self.tree_.n_node_samples[self.tree_.children_left[x]] +
                                       self.tree_.impurity[self.tree_.children_right[x]]*self.tree_.n_node_samples[self.tree_.children_right[x]]
                                      ) / (self.tree_.n_node_samples[x] - 2)
                f_score = f_score_numerator / f_score_denominator
                p = 1 - f.cdf(f_score, self.n_outputs_ * (self.tree_.n_node_samples[x] - 1), self.n_outputs_ * (self.tree_.n_node_samples[x] - 2))
                ps[x] = p
        self.p_value = ps

    def prune_node(self, index: int) -> None:
        """
        Prune all descendents of a given node. This node will become a leaf.

        Parameters
        ----------
        index
            Index of the node/leaf in the tree structure.

        Returns
        -------
        None
            Modified tree with all descendents of a given node pruned.
        """
        if self.is_node(index):
            stack = [index]
            while len(stack) > 0:
                 index_p = stack.pop()
                 if self.is_node(index_p):
                     stack.append(self.tree_.children_left[index_p])
                     stack.append(self.tree_.children_right[index_p])
                     self.tree_.children_left[index_p] = TREE_LEAF
                     self.tree_.children_right[index_p] = TREE_LEAF
                     self.tree_.feature[index_p] = TREE_UNDEFINED
                     self.tree_.threshold[index_p] = TREE_UNDEFINED

    def prune_tree(self, p_thres: float = 0.05) -> None:
        """
        Prune a tree based on F-test p values.

        Parameters
        ----------
        p_thres
            p-value threshold to prune nodes.
            (Default: `0.05`)

        Returns
        -------
        None
            Modified tree with unnecessary splits removed.
        """
        if p_thres < 0 or p_thres > 1:
            raise ValueError(
                    f"🛑 Please provide the `p_thres` between 0 and 1")
        stack = [0]
        while len(stack) > 0:
            index_n = stack.pop()
            if self.is_node(index_n):
                if self.p_value[index_n] < p_thres:
                    stack.append(self.tree_.children_left[index_n])
                    stack.append(self.tree_.children_right[index_n])
                else:
                    self.prune_node(index_n)

    def score(self, X, y, sample_weight = None) -> float:
        """
        Calculate the coefficient of determination between the prediction and truth.
        Different from multi-output problem where each output is calculated separately and the final R2 score is averaged across outputs, the score here is defined by considering each sample vector as a 'real' sample.

        Parameters
        ----------
        X
            Sample-by-feature query matrix.
        y
            Sample-by-output truth matrix.
        sample_weight
            Sample weights applied to squared distance of each sample.

        Returns
        ----------
        float
            Coefficient of determination.
        """
        if X.shape[1] != self.n_features_in_:
            raise ValueError(
                    f"🛑 Please provide `X` with {self.n_features_in_} columns")
        if not isinstance(y, np.ndarray):
            raise TypeError(
                    f"🛑 Please provide `y` as a numpy 2D array")
        if (X.shape[0] != y.shape[0]) or (self.n_outputs_ != y.shape[1]):
            raise ValueError(
                    f"🛑 Please provide `y` with {X.shape[0]} rows and {self.n_outputs_} columns")
        y_pred = self.predict(X)
        if sample_weight is not None:
            sample_weight = np.array(sample_weight)
            weight = sample_weight[:, np.newaxis]
        else:
            weight = 1.0
        numerator = (weight * (y - y_pred) ** 2).sum()
        denominator = (weight * (y - np.average(y, axis=0, weights = sample_weight)) ** 2).sum()
        if numerator == 0:
            return 1.0
        elif denominator == 0:
            return 0.0
        else:
            return 1 - numerator / denominator
