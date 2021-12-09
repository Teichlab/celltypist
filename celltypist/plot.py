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
    if use_as_prediction not in predictions.predicted_labels:
        if use_as_prediction == 'majority_voting':
            raise ValueError(f"🛑 Did not find the column `majority_voting` in the `AnnotationResult`, perform majority voting beforehand or use `use_as_prediction = 'predicted_labels'` instead")
        else:
            raise ValueError(f"🛑 Did not find such column '{use_as_prediction}', should be one of `majority_voting` or `predicted_labels`")
    pred = predictions.predicted_labels[use_as_prediction]
    #reference
    if isinstance(use_as_reference, str):
        if use_as_reference not in predictions.adata.obs:
            raise ValueError(f"🛑 Did not find such column '{use_as_reference}', please provide a valid metadata column")
        refer = predictions.adata.obs[use_as_reference]
    else:
        refer = np.array(use_as_reference)
        if len(refer) != len(pred):
            raise ValueError(f"🛑 Length of `use_as_reference` ({len(refer)}) provided does not match the number of cells ({len(pred)})")
    #score
    score = [row[pred[index]] for index, row in predictions.probability_matrix.iterrows()]
    #df x 2
    df = pd.DataFrame(dict(pred = pred, refer = refer, score = score))
    dot_size_df = df.pivot_table(values = 'score', index = 'pred', columns = 'refer', aggfunc = len, fill_value = 0, dropna = False, observed = True)
    dot_size_df = dot_size_df / dot_size_df.sum(axis = 0).values
    dot_color_df = df.pivot_table(values = 'score', index = 'pred', columns = 'refer', aggfunc = 'mean', fill_value = 0, dropna = False, observed = True)
    #reorder
    if prediction_order is None:
        prediction_order = pred.cat.categories
    else:
        if not np.all(np.unique(prediction_order) == np.unique(dot_size_df.index)):
            raise ValueError(f"🛑 Please provide a correct and comprehensive list of prediction cell types")
        prediction_order = np.array(prediction_order)
    dot_size_df = dot_size_df.loc[prediction_order]
    dot_color_df = dot_color_df.loc[prediction_order]
    if reference_order is None:
        reference_max_pred = dot_size_df.idxmax(axis = 0)
        reference_max_score = dot_size_df.max(axis = 0)
        sort_df = pd.DataFrame(dict(reference_order = dot_size_df.columns, reference_max_pred = reference_max_pred, reference_max_score = reference_max_score))
        sort_df['reference_max_pred'] = sort_df.reference_max_pred.astype('category')
        sort_df.reference_max_pred.cat.categories = [x for x in dot_size_df.index if x in sort_df.reference_max_pred.cat.categories]
        reference_order = sort_df.sort_values(by=['reference_max_pred', 'reference_max_score'], ascending = [True, False]).reference_order.values
    else:
        if not np.all(np.unique(reference_order) == np.unique(dot_size_df.columns)):
            raise ValueError(f"🛑 Please provide a correct and comprehensive list of reference cell types/clusters")
        reference_order = np.array(reference_order)
    dot_size_df = dot_size_df[reference_order]
    dot_color_df = dot_color_df[reference_order]
    #return
    return dot_size_df, dot_color_df
