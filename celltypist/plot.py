from .classifier import AnnotationResult
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
