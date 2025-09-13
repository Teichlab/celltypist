from .classifier import AnnotationResult
from .tree import TreeNode, Tree
import matplotlib
from matplotlib import pyplot as plt
from typing import Union, Optional
import numpy as np
import pandas as pd
import scanpy as sc

def _get_fraction_prob_df(predictions: AnnotationResult,
                          use_as_reference: Union[str, list, tuple, np.ndarray, pd.Series, pd.Index],
                          use_as_prediction: str = 'majority_voting',
                          prediction_order: Optional[Union[list, tuple, np.ndarray, pd.Series, pd.Index]] = None,
                          reference_order: Optional[Union[list, tuple, np.ndarray, pd.Series, pd.Index]] = None
                          ) -> tuple:
    """
    For internal use. Get the fraction and avg. probability data frames (predictions * truths) from AnnotationResult.
    """
    #prediction
    if not isinstance(predictions, AnnotationResult):
        raise TypeError(
                f"🛑 Please provide a correct input - an `AnnotationResult` derived from `celltypist.annotate`")
    if use_as_prediction not in predictions.predicted_labels:
        if use_as_prediction == 'majority_voting':
            raise KeyError(
                    f"🛑 Did not find the column `majority_voting` in the `AnnotationResult.predicted_labels`, perform majority voting beforehand or use `use_as_prediction = 'predicted_labels'` instead")
        else:
            raise KeyError(
                    f"🛑 Did not find such column '{use_as_prediction}', should be one of `'majority_voting'` or `'predicted_labels'`")
    pred = predictions.predicted_labels[use_as_prediction]
    #reference
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
    #score
    score = [(row[pred[index]] if pred[index] in row.index else row.max()) for index, row in predictions.probability_matrix.iterrows()]
    #df x 2
    df = pd.DataFrame(dict(pred = pred, refer = refer, score = score))
    dot_size_df = df.pivot_table(values = 'score', index = 'pred', columns = 'refer', aggfunc = len, fill_value = 0, dropna = False, observed = True)
    dot_size_df = dot_size_df / dot_size_df.sum(axis = 0).values
    dot_color_df = df.pivot_table(values = 'score', index = 'pred', columns = 'refer', aggfunc = 'mean', fill_value = 0, dropna = False, observed = True)
    #reorder
    if prediction_order is None:
        prediction_order = pred.cat.categories
    else:
        if not np.array_equal(np.sort(prediction_order), np.sort(dot_size_df.index)):
            raise ValueError(
                    f"🛑 Please provide a correct and comprehensive list of prediction cell types")
        prediction_order = np.array(prediction_order)
    dot_size_df = dot_size_df.loc[prediction_order]
    dot_color_df = dot_color_df.loc[prediction_order]
    if reference_order is None:
        reference_max_pred = dot_size_df.idxmax(axis = 0)
        reference_max_score = dot_size_df.max(axis = 0)
        sort_df = pd.DataFrame(dict(reference_order = dot_size_df.columns, reference_max_pred = reference_max_pred, reference_max_score = reference_max_score))
        sort_df['reference_max_pred'] = sort_df.reference_max_pred.astype('category')
        sort_df['reference_max_pred'] = sort_df.reference_max_pred.cat.reorder_categories([x for x in dot_size_df.index if x in sort_df.reference_max_pred.cat.categories])
        reference_order = sort_df.sort_values(by=['reference_max_pred', 'reference_max_score'], ascending = [True, False]).reference_order.values
    else:
        if not np.array_equal(np.sort(reference_order), np.sort(dot_size_df.columns)):
            raise ValueError(
                    f"🛑 Please provide a correct and comprehensive list of reference cell types/clusters")
        reference_order = np.array(reference_order)
    dot_size_df = dot_size_df[reference_order]
    dot_color_df = dot_color_df[reference_order]
    #return
    return dot_size_df, dot_color_df

def dotplot(
            #get size and color df
            predictions: AnnotationResult,
            use_as_reference: Union[str, list, tuple, np.ndarray, pd.Series, pd.Index],
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
    use_as_reference
        Key (column name) of the input AnnData representing the reference cell types (or clusters) celltypist will assess.
        Also accepts any list-like objects already loaded in memory (such as an array).
    use_as_prediction
        Column name of :attr:`~celltypist.classifier.AnnotationResult.predicted_labels` specifying the prediction type which the assessment is based on.
        Set to `'predicted_labels'` if you want to assess the prediction result without majority voting.
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
    #df x 2
    dot_size_df, dot_color_df = _get_fraction_prob_df(predictions, use_as_reference, use_as_prediction, None, None)
    #reference
    reference_order = reference_order if reference_order is not None else dot_size_df.columns
    #prediction
    if prediction_order is None:
        if filter_prediction < 0 or filter_prediction > 1:
            raise ValueError(
                    f"🛑 Please provide the `filter_prediction` between 0 and 1")
        keep_pred = dot_size_df.max(axis = 1) >= filter_prediction
        prediction_order = dot_size_df.index[keep_pred]
    #in case reference_order or prediction_order is string
    if isinstance(reference_order, str):
        reference_order = [reference_order]
    if isinstance(prediction_order, str):
        prediction_order = [prediction_order]
    #subset
    dot_size_df = dot_size_df.loc[prediction_order, reference_order]
    dot_color_df = dot_color_df.loc[prediction_order, reference_order]
    #column to string
    dot_size_df.columns = dot_size_df.columns.astype(str)
    dot_color_df.columns = dot_color_df.columns.astype(str)
    #AnnData, groupby, and var_names
    _adata = sc.AnnData(np.zeros(dot_size_df.shape))
    _adata.var_names = dot_size_df.columns
    _adata.obs_names = dot_size_df.index
    _adata.obs['_pred'] = dot_size_df.index
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

def _assign_coords(node: TreeNode, depth: int = 1, y: Optional[list] = None, coords: Optional[dict] = None, type: str = "cladogram") -> dict:
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
            _assign_coords(child, depth + 1, y, coords, type)
        if type == "phylogram":
            child_ordinates = [coords[child.original_name][1] for child in node.children]
            ordinate = sum(child_ordinates) / len(child_ordinates)
        else:
            leaf_ordinates = [coords[leaf][1] for leaf in node.cell_types(leaf_only = True)]
            ordinate = sum(leaf_ordinates) / len(leaf_ordinates)
        coords[node.original_name] = (depth - 1, ordinate, False)
    return coords

def treeviz(tree: Tree,
        #design
        type: str = "cladogram", direction: str = "right", sort: bool = False, recursive: bool = True, descending: bool = True,
        #branch
        edge_color: str = '#0000007B', edge_width: Optional[float] = None,
        #node
        node_shape: str = "o", node_color: str = '#2E91E5', node_size: Optional[float] = None,
        #leaf
        leaf_shape: Optional[str] = None, leaf_color: Optional[str] = None, leaf_size: Optional[float] = None,
        #show
        show_node_label: bool = False, show_leaf_label: bool = True,
        #node label
        node_label_color: str = '#000000', node_label_size: Optional[Union[float, str]] = None, node_label_ha: str = "center", node_label_va: str = "bottom", node_label_rotation: Optional[Union[float, str]] = None,
        #leaf label
        leaf_label_color: Optional[str] = None, leaf_label_size: Optional[Union[float, str]] = None, leaf_label_ha: Optional[str] = None, leaf_label_va: Optional[str] = None, leaf_label_rotation: Optional[Union[float, str]] = None,
        #figure elements
        title: Optional[str] = None,
        #show and/or save figure
        ax: Optional[matplotlib.axes.Axes] = None, figsize: Optional[Union[list, tuple]] = None, show: bool = True, save: Union[str, bool] = False,
        #others
        edge_dict = {}, node_dict = {}, leaf_dict = {}, node_label_dict = {}, leaf_label_dict = {},
        ) -> None:
    """
    Visualize a cell type hierarchical tree.

    Parameters
    ----------
    tree
        A :class:`~celltypist.tree.Tree` object to visualize.
    type
        Layout style of the tree, either `'cladogram'` or `'phylogram'`.
        (Default: `'cladogram'`)
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
    if type not in ("cladogram", "phylogram"):
        raise ValueError(
                f"🛑 `type` must be 'cladogram' or 'phylogram'")
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
    leaf_label_color =  node_label_color if leaf_label_color is None else leaf_label_color
    leaf_label_size = node_label_size if leaf_label_size is None else leaf_label_size
    #coords
    if sort:
        tree = tree.copy()
        tree.sort_tree(recursive = recursive, descending = descending)
    coords = _assign_coords(tree.root, type = type)
    oriented_coords = {}
    if direction == "right":
        for name, (x, y, is_leaf) in coords.items():
            oriented_coords[name] = (x, y, is_leaf)
    else:
        max_new_y = tree.depth - 1
        for name, (x, y, is_leaf) in coords.items():
            new_y = max_new_y - x
            oriented_coords[name] = (y, new_y, is_leaf)
    #axes
    if ax is None:
        if figsize is None:
            figsize = (tree.depth * 2.5, tree.n_leaves * 0.6) if direction == "right" else (tree.n_leaves * 0.6, tree.depth * 2.5)
        _, ax = plt.subplots(figsize = figsize)
    #edges
    if type == "cladogram":
        for node, (x, y, _) in oriented_coords.items():
            parent = tree.find_parent(node)
            if parent is not None:
                xp, yp, _ = oriented_coords[parent.original_name]
                ax.plot([xp, x], [yp, y], color = edge_color, lw = edge_width, **edge_dict)
    else:
        def draw_phylo(node):
            if node.is_leaf():
                return
            x, y, _ = oriented_coords[node.original_name]
            child_coords = [oriented_coords[c.original_name] for c in node.children]
            child_xs = [cx for cx, _, _ in child_coords]
            child_ys = [cy for _, cy, _ in child_coords]
            if direction == "right":
                ax.plot([x, x], [min(child_ys), max(child_ys)], color = edge_color, lw = edge_width, **edge_dict)
                for (cx, cy, _) in child_coords:
                    ax.plot([x, cx], [cy, cy], color = edge_color, lw = edge_width, **edge_dict)
            else:
                ax.plot([min(child_xs), max(child_xs)], [y, y], color = edge_color, lw = edge_width, **edge_dict)
                for (cx, cy, _) in child_coords:
                    ax.plot([cx, cx], [y, cy], color = edge_color, lw = edge_width, **edge_dict)
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
        ax.set(xlim = [-0.5, tree.depth], ylim = [-0.5, tree.n_leaves - 0.5], title = title)
    else:
        ax.set(xlim = [-0.5, tree.n_leaves - 0.5], ylim = [-1, tree.depth - 0.5], title = title)
    ax.set_axis_off()
    #show and save
    if save:
        plt.savefig(save) if isinstance(save, str) else plt.savefig('treeviz.pdf')
    if show:
        plt.show()
    if save:
        plt.close()
