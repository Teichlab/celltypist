from .classifier import AnnotationResult
from .tree import TreeNode, Tree
import matplotlib
from matplotlib import pyplot as plt
from typing import Union, Optional
import numpy as np
import pandas as pd
import scanpy as sc

def _get_fraction_prob_df(df: pd.DataFrame,
                          prediction_order: Optional[Union[str, list, tuple, np.ndarray, pd.Series, pd.Index]] = None,
                          reference_order: Optional[Union[str, list, tuple, np.ndarray, pd.Series, pd.Index]] = None,
                          filter_prediction: float = 0.0) -> tuple:
    """
    For internal use. Get the fraction and avg. probability data frames (predictions * truths) from the prediction-truth-score data frame.
    """
    if not all(isinstance(x, (int, float)) for x in df.iloc[:, 2]):
        raise TypeError(
                f"🛑 The third column '{df.columns[2]}' must contain only floats or integers")
    df = df.copy()
    df.columns = ['pred', 'refer', 'score']
    prediction_order_all = None
    if isinstance(df.pred.dtype, pd.CategoricalDtype):
        prediction_order_all = df.pred.cat.remove_unused_categories().cat.categories
    for col in ('pred', 'refer'):
        df[col] = df[col].astype(str)
    #df x 2
    dot_size_df = df.pivot_table(values = 'score', index = 'pred', columns = 'refer', aggfunc = len, fill_value = 0, dropna = False, observed = True)
    dot_size_df = dot_size_df / dot_size_df.sum(axis = 0).values
    dot_color_df = df.pivot_table(values = 'score', index = 'pred', columns = 'refer', aggfunc = 'mean', fill_value = 0, dropna = False, observed = True)
    #reorder
    if prediction_order_all is not None:
        dot_size_df = dot_size_df.loc[prediction_order_all]
        dot_color_df = dot_color_df.loc[prediction_order_all]
    reference_max_pred = dot_size_df.idxmax(axis = 0)
    reference_max_score = dot_size_df.max(axis = 0)
    sort_df = pd.DataFrame(dict(reference_order = dot_size_df.columns, reference_max_pred = reference_max_pred, reference_max_score = reference_max_score))
    sort_df['reference_max_pred'] = sort_df.reference_max_pred.astype('category')
    sort_df['reference_max_pred'] = sort_df.reference_max_pred.cat.reorder_categories([x for x in dot_size_df.index if x in sort_df.reference_max_pred.cat.categories])
    reference_order_all = sort_df.sort_values(by=['reference_max_pred', 'reference_max_score'], ascending = [True, False]).reference_order.values
    dot_size_df = dot_size_df[reference_order_all]
    dot_color_df = dot_color_df[reference_order_all]
    #subset of order
    reference_order = dot_size_df.columns if reference_order is None else reference_order
    if isinstance(reference_order, str):
        reference_order = [reference_order]
    if prediction_order is None:
        if filter_prediction < 0 or filter_prediction > 1:
            raise ValueError(
                    f"🛑 Please provide the `filter_prediction` between 0 and 1")
        keep_pred = dot_size_df.max(axis = 1) >= filter_prediction
        prediction_order = dot_size_df.index[keep_pred]
    if isinstance(prediction_order, str):
        prediction_order = [prediction_order]
    dot_size_df = dot_size_df.loc[prediction_order, reference_order]
    dot_color_df = dot_color_df.loc[prediction_order, reference_order]
    #return
    return dot_size_df, dot_color_df

def dotplot(
            #get size and color df
            predictions: Union[AnnotationResult, pd.DataFrame, tuple, list],
            use_as_reference: Optional[Union[str, list, tuple, np.ndarray, pd.Series, pd.Index]] = None,
            use_as_prediction: str = 'majority_voting',
            prediction_order: Optional[Union[str, list, tuple, np.ndarray, pd.Series, pd.Index]] = None,
            reference_order: Optional[Union[str, list, tuple, np.ndarray, pd.Series, pd.Index]] = None,
            filter_prediction: float = 0.0,
            #color
            cmap: str = 'RdBu_r',
            vmin: Optional[float] = 0.0,
            vmax: Optional[float] = 1.0,
            #color bar
            colorbar_title: Optional[str] = 'Mean probability',
            #size
            dot_min: Optional[float] = 0.0,
            dot_max: Optional[float] = 1.0,
            smallest_dot: Optional[float] = 0.0,
            #size bar
            size_title: Optional[str] = 'Fraction of cells (%)',
            #global
            swap_axes: Optional[bool] = False,
            title: Optional[str] = 'CellTypist label transfer',
            figsize: Optional[tuple] = None,
            #display
            show: Optional[bool] = None,
            save: Union[str, bool, None] = None,
            ax: Optional[sc.pl._utils._AxesSubplot] = None,
            return_fig: Optional[bool] = False,
            #other
            **kwds
           ) -> Union[sc.pl.DotPlot, dict, None]:
    """
    Generate a dot plot showing CellTypist label transfer. This is a wrapper around the :class:`scanpy.pl.DotPlot` with selected parameters and customized defaults.

    Parameters
    ----------
    predictions
        An :class:`~celltypist.classifier.AnnotationResult` object containing celltypist prediction result through :func:`~celltypist.annotate`.
        Can also be a :class:`~pandas.DataFrame` whose columns are ordered as prediction, truth, and score, or a tuple/list containing (dot_size_df, dot_color_df).
    use_as_reference
        Key (column name) of the input AnnData representing the reference cell types (or clusters) celltypist will assess.
        Also accepts any list-like objects already loaded in memory (such as an array).
        This argument is valid only when the input is an :class:`~celltypist.classifier.AnnotationResult` object.
    use_as_prediction
        Column name of :attr:`~celltypist.classifier.AnnotationResult.predicted_labels` specifying the prediction type which the assessment is based on.
        Set to `'predicted_labels'` if you want to assess the prediction result without majority voting.
        This argument is valid only when the input is an :class:`~celltypist.classifier.AnnotationResult` object.
        (Default: `'majority_voting'`)
    prediction_order
        Order in which to show the predicted cell types. Can be a subset of predicted cell type labels.
        Default to plotting all predicted labels, with the order of categories as is (alphabetical order in most cases).
    reference_order
        Order in which to show the reference cell types (or clusters). Can be a subset of reference cell types (or clusters).
        Default to plotting all reference cell types, with an order that ensures the resulting dot plot is diagonal.
    filter_prediction
        Filter out the predicted cell types with the maximal assignment fractions less than `filter_prediction`.
        This argument is only effective when `prediction_order` is not specified, and can be used to reduce the number of predicted cell types displayed in the dot plot.
        Default to 0 (no filtering).
    title
        Title of the dot plot.
        (Default: `'CellTypist label transfer'`)
    size_title
        Legend title for the dot sizes.
        (Default: `'Fraction of cells (%)'`)
    colorbar_title
        Legend title for the dot colors.
        (Default: `'Mean probability'`)
    swap_axes
        Whether to swap the x and y axes.
        (Default: `False`)
    others
        All other parameters are the same as :func:`scanpy.pl.dotplot` with selected tags and customized defaults.

    Returns
    ----------
    If `return_fig` is `True`, returns a :class:`scanpy.pl.DotPlot` object, else if `show` is false, return axes dict.
    """
    if isinstance(predictions, AnnotationResult):
        if use_as_prediction not in predictions.predicted_labels:
            if use_as_prediction == 'majority_voting':
                raise KeyError(
                        f"🛑 Did not find the column `majority_voting` in the `AnnotationResult.predicted_labels`, perform majority voting beforehand or use `use_as_prediction = 'predicted_labels'` instead")
            else:
                raise KeyError(
                        f"🛑 Did not find such column '{use_as_prediction}', should be one of `'majority_voting'` or `'predicted_labels'`")
        pred = predictions.predicted_labels[use_as_prediction]
        if use_as_reference is None:
            raise ValueError(
                    f"🛑 `use_as_reference` must be provided")
        if isinstance(use_as_reference, str):
            if use_as_reference not in predictions.adata.obs:
                raise KeyError(
                        f"🛑 Did not find such column '{use_as_reference}', please provide a valid metadata column")
            refer = predictions.adata.obs[use_as_reference]
        else:
            refer = np.array(use_as_reference)
            if len(refer) != len(pred):
                raise ValueError(
                        f"🛑 Length of `use_as_reference` ({len(refer)}) provided does not match the number of cells ({len(pred)})")
        score = [(row[pred[index]] if pred[index] in row.index else row.max()) for index, row in predictions.probability_matrix.iterrows()]
        df = pd.DataFrame(dict(pred = pred, refer = refer, score = score))
        dot_size_df, dot_color_df = _get_fraction_prob_df(df, prediction_order, reference_order, filter_prediction)
    elif isinstance(predictions, pd.DataFrame) and predictions.shape[1] == 3:
        dot_size_df, dot_color_df = _get_fraction_prob_df(predictions, prediction_order, reference_order, filter_prediction)
    elif isinstance(predictions, (tuple, list)) and len(predictions) == 2 and isinstance(predictions[0], pd.DataFrame) and isinstance(predictions[1], pd.DataFrame):
        dot_size_df, dot_color_df = predictions
    else:
        raise Exception(
                f"🛑 Invalid input")
    #AnnData, groupby, and var_names
    _adata = sc.AnnData(np.zeros(dot_size_df.shape))
    _adata.var_names = dot_size_df.columns
    _adata.obs_names = dot_size_df.index
    _adata.obs['_pred'] = dot_size_df.index
    _adata.obs['_pred'] = _adata.obs['_pred'].astype('category').cat.reorder_categories(dot_size_df.index)
    #DotPlot
    dp = sc.pl.DotPlot(_adata, dot_size_df.columns, '_pred', title = title, figsize = figsize, dot_color_df = dot_color_df, dot_size_df = dot_size_df, ax = ax, vmin = vmin, vmax = vmax, **kwds)
    if swap_axes:
        dp.swap_axes()
    dp = dp.style(cmap = cmap, dot_max = dot_max, dot_min = dot_min, smallest_dot = smallest_dot, dot_edge_lw = kwds.pop('linewidth', 0.2)).legend(colorbar_title = colorbar_title, size_title = size_title)
    if return_fig:
        return dp
    else:
        dp.make_figure()
        sc.pl._utils.savefig_or_show('CellTypist_dotplot_', show = show, save = save)
        show = sc._settings.settings.autoshow if show is None else show
        if not show:
            return dp.get_axes()

def _assign_coords(node: TreeNode, depth: int = 1, y: Optional[list] = None, coords: Optional[dict] = None, layout: str = "diagonal", extend_leaves: bool = False) -> dict:
    """
    For internal use. Assign coordinates (x, y, is_leaf) to a node and its descendants.
    """
    if coords is None:
        coords, y = {}, [0]
    if node.is_leaf():
        coords[node.original_name] = (depth - 1, y[0], True)
        y[0] += 1
    else:
        for child in node.children:
            _assign_coords(child, depth + 1, y, coords, layout, extend_leaves)
        if layout == "rectangular":
            child_ordinates = [coords[child.original_name][1] for child in node.children]
            ordinate = sum(child_ordinates) / len(child_ordinates)
        else:
            leaf_ordinates = [coords[leaf][1] for leaf in node.cell_types(leaf_only = True)]
            ordinate = sum(leaf_ordinates) / len(leaf_ordinates)
        coords[node.original_name] = (depth - 1, ordinate, False)
    if depth == 1 and extend_leaves:
        max_x = max(x for x, _, _ in coords.values())
        for name, (x, y, is_leaf) in coords.items():
            if is_leaf:
                coords[name] = (max_x, y, True)
    return coords

def treeviz(tree: Tree,
        #design
        layout: str = "diagonal", direction: str = "right", sort: bool = False, recursive: bool = True, descending: bool = True,
        #branch
        edge_color: str = '#0000007B', edge_width: Optional[float] = None,
        #node
        node_shape: str = "o", node_color: str = '#2E91E5', node_size: Optional[float] = None,
        #leaf
        leaf_shape: Optional[str] = None, leaf_color: Optional[str] = None, leaf_size: Optional[float] = None,
        #node and leaf color map
        node_color_map: Optional[dict] = None, cmap: Union[matplotlib.colors.Colormap, str] = 'Reds', cmap_min: Optional[float] = None, cmap_max: Optional[float] = None,
        #node and leaf size map
        node_size_map: Optional[dict] = None, smap: Optional[tuple] = None, smap_min: Optional[float] = None, smap_max: Optional[float] = None,
        #show
        show_node_label: bool = False, show_leaf_label: bool = True,
        #node label
        node_label_color: str = '#000000', node_label_size: Optional[Union[float, str]] = None, node_label_ha: str = "center", node_label_va: str = "bottom", node_label_rotation: Optional[Union[float, str]] = None,
        #leaf label
        leaf_label_color: Optional[str] = None, leaf_label_size: Optional[Union[float, str]] = None, leaf_label_ha: Optional[str] = None, leaf_label_va: Optional[str] = None, leaf_label_rotation: Optional[Union[float, str]] = None,
        #extend leaves
        extend_leaves: bool = False,
        #figure elements
        title: Optional[str] = None,
        #show and/or save figure
        ax: Optional[matplotlib.axes.Axes] = None, figsize: Optional[Union[list, tuple]] = None, show: bool = True, save: Union[str, bool] = False,
        #others
        edge_dict: Optional[dict] = None, node_dict: Optional[dict] = None, leaf_dict: Optional[dict] = None, node_label_dict: Optional[dict] = None, leaf_label_dict: Optional[dict] = None,
        ) -> None:
    """
    Visualize a cell type hierarchical tree.

    Parameters
    ----------
    tree
        A :class:`~celltypist.tree.Tree` object to visualize.
    layout
        Layout style of the tree, either `'diagonal'` or `'rectangular'`.
        (Default: `'diagonal'`)
    direction
        Direction of tree growth, either `'right'` or `'down'`.
        (Default: `'right'`)
    sort
        Whether to sort children using :meth:`~celltypist.tree.Tree.sort_tree` before plotting.
        (Default: `False`)
    recursive
        If `sort = True`, whether to sort by the total number of leaves recursively contained.
        (Default: `True`)
    descending
        If `sort = True`, whether to sort in descending order.
        (Default: `True`)
    edge_color
        Color of edges/branches.
        (Default: `'#0000007B'`)
    edge_width
        Width of edges/branches in points.
        Default to 1.5 in a canonical Matplotlib setting.
    node_shape
        Shape of internal nodes.
        (Default: `'o'`)
    node_color
        Color of internal nodes.
        (Default: `'#2E91E5'`)
    node_size
        Size of internal nodes in points.
        Default to 6.0 in a canonical Matplotlib setting.
    leaf_shape
        Shape of leaf nodes. Default to `node_shape`.
    leaf_color
        Color of leaf nodes. Default to `node_color`.
    leaf_size
        Size of leaf nodes in points. Default to `node_size`.
    node_color_map
        Optional mapping from node names (internal or leaf) to color specifications. Keys must be a subset of all node names in the tree.
        Values can be either valid color specifications which are used directly, or numeric values which are normalized by `cmap_min` and `cmap_max` and then mapped to `cmap`.
        Nodes not present in `node_color_map` fall back to `node_color` for internal nodes and `leaf_color` for leaf nodes.
    cmap
        Colormap used when `node_color_map` contains continuous numeric values.
        (Default: `'Reds'`)
    cmap_min
        Lower bound for color normalization when `node_color_map` contains numeric values. Values outside the normalization range are clipped.
        Default to the minimum of the provided values.
    cmap_max
        Upper bound for color normalization when `node_color_map` contains numeric values. Values outside the normalization range are clipped.
        Default to the maximum of the provided values.
    node_size_map
        Optional mapping from node names (internal or leaf) to continuous numeric values that control marker size. Keys must be a subset of all node names in the tree.
        Values are normalized by `smap_min` and `smap_max` and then mapped to `smap`.
        Nodes not present in `node_size_map` fall back to `node_size` for internal nodes and `leaf_size` for leaf nodes.
    smap
        Two-element sequence specifying the minimum and maximum marker sizes.
        Default to (default_marker_size / 2, default_marker_size * 2), which is (3, 12) in a canonical Matplotlib setting.
    smap_min
        Lower bound for size normalization when `node_size_map` is provided. Values outside the normalization range are clipped.
        Default to the minimum of the provided values.
    smap_max
        Upper bound for size normalization when `node_size_map` is provided. Values outside the normalization range are clipped.
        Default to the maximum of the provided values.
    show_node_label
        Whether to show labels for internal nodes.
        (Default: `False`)
    show_leaf_label
        Whether to show labels for leaf nodes.
        (Default: `True`)
    node_label_color
        Color of internal node labels.
        (Default: `'#000000'`)
    node_label_size
        Size of internal node labels.
        Default to 10.0 in a canonical Matplotlib setting.
    node_label_ha
        Horizontal alignment of internal node labels.
        (Default: `'center'`)
    node_label_va
        Vertical alignment of internal node labels.
        (Default: `'bottom'`)
    node_label_rotation
        Rotation angle of internal node labels.
        Default to 0.0 (no rotation) in a canonical Matplotlib setting.
    leaf_label_color
        Color of leaf labels. Default to `node_label_color`.
    leaf_label_size
        Size of leaf labels. Default to `node_label_size`.
    leaf_label_ha
        Horizontal alignment of leaf labels. Auto-set by `direction` if not provided.
    leaf_label_va
        Vertical alignment of leaf labels. Auto-set by `direction` if not provided.
    leaf_label_rotation
        Rotation angle of leaf labels. Auto-set by `direction` if not provided.
    extend_leaves
        Whether to horizontally extend all leaf nodes to the maximal depth of the tree.
        (Default: `False`)
    title
        Figure title. Default to `"Cell type tree: {tree.handle}"`.
    ax
        An :class:`~matplotlib.axes.Axes` where the tree will be drawn.
        Default to draw on a new axes.
    figsize
        Tuple of figure width and height in inches.
        Default to auto-adjust based on tree depth and number of leaves.
    show
        Whether to display the figure.
        (Default: `True`)
    save
        Whether to save the figure. This can also be a figure filename.
        (Default: `False`)
    edge_dict
        Extra keyword arguments passed to :class:`~matplotlib.lines.Line2D` for edges/branches.
    node_dict
        Extra keyword arguments passed to :class:`~matplotlib.lines.Line2D` for internal nodes.
    leaf_dict
        Extra keyword arguments passed to :class:`~matplotlib.lines.Line2D` for leaf nodes.
    node_label_dict
        Extra keyword arguments passed to :class:`~matplotlib.text.Text` for internal node labels.
    leaf_label_dict
        Extra keyword arguments passed to :class:`~matplotlib.text.Text` for leaf labels.

    Returns
    ----------
    None
    """
    #params
    if not isinstance(tree, Tree):
        raise TypeError(
                f"🛑 Please provide a `Tree` instance")
    if layout not in ("diagonal", "rectangular"):
        raise ValueError(
                f"🛑 `layout` must be 'diagonal' or 'rectangular'")
    if direction == "right":
        leaf_label_ha = "left" if leaf_label_ha is None else leaf_label_ha
        leaf_label_va = "center" if leaf_label_va is None else leaf_label_va
        leaf_label_rotation = 0 if leaf_label_rotation is None else leaf_label_rotation
    elif direction == "down":
        leaf_label_ha = "center" if leaf_label_ha is None else leaf_label_ha
        leaf_label_va = "top" if leaf_label_va is None else leaf_label_va
        leaf_label_rotation = 90 if leaf_label_rotation is None else leaf_label_rotation
    else:
        raise ValueError(
                f"🛑 `direction` must be 'right' or 'down'")
    leaf_shape = node_shape if leaf_shape is None else leaf_shape
    leaf_color = node_color if leaf_color is None else leaf_color
    leaf_size = node_size if leaf_size is None else leaf_size
    leaf_label_color = node_label_color if leaf_label_color is None else leaf_label_color
    leaf_label_size = node_label_size if leaf_label_size is None else leaf_label_size
    edge_dict = edge_dict or {}
    node_dict = node_dict or {}
    leaf_dict = leaf_dict or {}
    node_label_dict = node_label_dict or {}
    leaf_label_dict = leaf_label_dict or {}
    #coords
    if sort:
        tree = tree.copy()
        tree.sort_tree(recursive = recursive, descending = descending)
    coords = _assign_coords(tree.root, layout = layout, extend_leaves = extend_leaves)
    tree_depth = tree.depth
    oriented_coords = {}
    if direction == "right":
        for name, (x, y, is_leaf) in coords.items():
            oriented_coords[name] = (x, y, is_leaf)
    else:
        for name, (x, y, is_leaf) in coords.items():
            oriented_coords[name] = (y, tree_depth - 1 - x, is_leaf)
    #axes
    tree_n_leaves = tree.n_leaves
    if ax is None:
        if figsize is None:
            figsize = (tree_depth * 2.5, tree_n_leaves * 0.6) if direction == "right" else (tree_n_leaves * 0.6, tree_depth * 2.5)
        _, ax = plt.subplots(figsize = figsize)
    #edges
    if layout == "diagonal":
        for node, (x, y, _) in oriented_coords.items():
            parent = tree.find_parent(node)
            if parent is not None:
                xp, yp, _ = oriented_coords[parent.original_name]
                ax.plot([xp, x], [yp, y], color = edge_color, lw = edge_width, marker = 'None', **edge_dict)
    else:
        def draw_phylo(node):
            if node.is_leaf():
                return
            x, y, _ = oriented_coords[node.original_name]
            child_coords = [oriented_coords[c.original_name] for c in node.children]
            child_xs = [cx for cx, _, _ in child_coords]
            child_ys = [cy for _, cy, _ in child_coords]
            if direction == "right":
                ax.plot([x, x], [min(child_ys), max(child_ys)], color = edge_color, lw = edge_width, marker = 'None', **edge_dict)
                for (cx, cy, _) in child_coords:
                    ax.plot([x, cx], [cy, cy], color = edge_color, lw = edge_width, marker = 'None', **edge_dict)
            else:
                ax.plot([min(child_xs), max(child_xs)], [y, y], color = edge_color, lw = edge_width, marker = 'None', **edge_dict)
                for (cx, cy, _) in child_coords:
                    ax.plot([cx, cx], [y, cy], color = edge_color, lw = edge_width, marker = 'None', **edge_dict)
            for c in node.children:
                draw_phylo(c)
        draw_phylo(tree.root)
    #nodes & labels
    for name, (x, y, is_leaf) in oriented_coords.items():
        if is_leaf:
            ax.plot(x, y, marker = leaf_shape, ms = leaf_size, color = leaf_color, ls = 'None', **leaf_dict)
            if show_leaf_label:
                ax.text(x + 0.05 if direction == "right" else x, y - 0.05 if direction == "down" else y, name, color = leaf_label_color, size = leaf_label_size, ha = leaf_label_ha, va = leaf_label_va, rotation = leaf_label_rotation, **leaf_label_dict)
        else:
            ax.plot(x, y, marker = node_shape, ms = node_size, color = node_color, ls = 'None', **node_dict)
            if show_node_label:
                ax.text(x, y + 0.05, name, color = node_label_color, size = node_label_size, ha = node_label_ha, va = node_label_va, rotation = node_label_rotation, **node_label_dict)
    #frame
    title = f"Cell type tree: {tree.handle}" if title is None else title
    if direction == "right":
        ax.set(xlim = [-0.5, tree_depth], ylim = [-0.5, tree_n_leaves - 0.5], title = title)
    else:
        ax.set(xlim = [-0.5, tree_n_leaves - 0.5], ylim = [-1, tree_depth - 0.5], title = title)
    ax.set_axis_off()
    #show and save
    if save:
        plt.savefig(save) if isinstance(save, str) else plt.savefig('CellTypist_treeviz.pdf')
    if show:
        plt.show()
    if save:
        plt.close()

treevis = treeviz
