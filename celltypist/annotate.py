from . import classifier, models
from typing import Optional, Union
import numpy as np
import pandas as pd

def annotate(filename: str,
             model: str = "",
             transpose_input: bool = False,
             gene_file: Optional[str] = None,
             cell_file: Optional[str] = None,
             majority_voting: bool = False,
             over_clustering: Optional[Union[str, list, np.ndarray, pd.Series]] = None) -> classifier.AnnotationResult:
    """
    Run the prediction and (optional) majority voting to annotate the input dataset.

    Parameters
    ----------
    filename
        Path to the input count matrix (supported types are csv, txt, tsv, tab and mtx) or Scanpy object (h5ad).
        If it's the former, a cell-by-gene format is desirable (see `transpose_input` for more information).
        Genes should be gene symbols. Non-expressed genes are preferred to be provided as well.
    model
        Model used to predict the input cells. Default to using the `Immune_All_Low.pkl` model.
        To see all available models and their descriptions, use :func:`~celltypist.models.models_description()`.
    transpose_input
        Whether to transpose the input matrix. Set to `True` if `filename` is provided in a gene-by-cell format.
        (Default: `False`)
    gene_file
        Path to the file which stores each gene per line corresponding to the genes used in the provided mtx file.
        Ignored if `filename` is not provided in the mtx format.
    cell_file
        Path to the file which stores each cell per line corresponding to the cells used in the provided mtx file.
        Ignored if `filename` is not provided in the mtx format.
    majority_voting
        Whether to refine the predicted labels by running the majority voting classifier after over-clustering.
        (Default: `False`)
    over_clustering
        This argument can be provided in several ways:
        1) an input plain file with the over-clustering result of one cell per line.
        2) a string key specifying an existing metadata column in the `AnnData` (pre-created by the user).
        3) a python list, numpy array, or pandas series indicating the over-clustering result of all cells.
        4) if none of the above is provided, will use a heuristic over-clustering approach based on input data size.
        Ignored if `majority_voting` is set to `False`.

    Returns
    ----------
    :class:`~celltypist.classifier.AnnotationResult`
        An :class:`~celltypist.classifier.AnnotationResult` object. Two important attributes within this class are:
        1) :attr:`~celltypist.classifier.AnnotationResult.predicted_labels`, predicted labels from celltypist.
        2) :attr:`~celltypist.classifier.AnnotationResult.probability_table`, probability matrix from celltypist.
    """
    #load model
    sgd_classifier = models.load(model)
    #construct Classifier class
    clf = classifier.Classifier(filename = filename, model = sgd_classifier, transpose = transpose_input, gene_file = gene_file, cell_file = cell_file)
    #predict
    predictions = clf.celltype()
    if majority_voting is False:
        return predictions
    #over clustering
    if over_clustering is None:
        over_clustering = clf.over_cluster()
    elif isinstance(over_clustering, str):
        if over_clustering in clf.adata.obs:
            over_clustering = clf.adata.obs[over_clustering]
        else:
            try:
                with open(over_clustering, 'rt') as f:
                    over_clustering = [x.strip() for x in f.readlines()]
            except Exception as e:
                raise Exception(f"🛑 {e}")
    if len(over_clustering) != clf.adata.shape[0]:
        raise ValueError(f"🛑 Length of `over_clustering` ({len(over_clustering)}) does not match the number of input cells ({clf.adata.shape[0]})")
    #majority voting
    return classifier.Classifier.majority_vote(predictions, over_clustering)
