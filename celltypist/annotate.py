from . import classifier
from .models import Model
from typing import Optional, Union
import numpy as np
import pandas as pd
from anndata import AnnData
from . import logger

def annotate(filename: Union[AnnData,str] = "",
             model: Optional[Union[str, Model]] = None,
             transpose_input: bool = False,
             gene_file: Optional[str] = None,
             cell_file: Optional[str] = None,
             mode: str = 'best match',
             p_thres: float = 0.5,
             majority_voting: bool = False,
             over_clustering: Optional[Union[str, list, tuple, np.ndarray, pd.Series, pd.Index]] = None,
             use_GPU: bool = False,
             min_prop: float = 0) -> classifier.AnnotationResult:
    """
    Run the prediction and (optional) majority voting to annotate the input dataset.

    Parameters
    ----------
    filename
        Path to the input count matrix (supported types are csv, txt, tsv, tab and mtx) or AnnData (h5ad).
        If it's the former, a cell-by-gene format is desirable (see `transpose_input` for more information).
        Also accepts the input as an :class:`~anndata.AnnData` object already loaded in memory.
        Genes should be gene symbols. Non-expressed genes are preferred to be provided as well.
    model
        Model used to predict the input cells. Default to using the `'Immune_All_Low.pkl'` model.
        Can be a :class:`~celltypist.models.Model` object that wraps the logistic Classifier and the StandardScaler, the
        path to the desired model file, or the model name.
        To see all available models and their descriptions, use :func:`~celltypist.models.models_description`.
    transpose_input
        Whether to transpose the input matrix. Set to `True` if `filename` is provided in a gene-by-cell format.
        (Default: `False`)
    gene_file
        Path to the file which stores each gene per line corresponding to the genes used in the provided mtx file.
        Ignored if `filename` is not provided in the mtx format.
    cell_file
        Path to the file which stores each cell per line corresponding to the cells used in the provided mtx file.
        Ignored if `filename` is not provided in the mtx format.
    mode
        The way cell prediction is performed.
        For each query cell, the default (`'best match'`) is to choose the cell type with the largest score/probability as the final prediction.
        Setting to `'prob match'` will enable a multi-label classification, which assigns 0 (i.e., unassigned), 1, or >=2 cell type labels to each query cell.
        (Default: `'best match'`)
    p_thres
        Probability threshold for the multi-label classification. Ignored if `mode` is `'best match'`.
        (Default: 0.5)
    majority_voting
        Whether to refine the predicted labels by running the majority voting classifier after over-clustering.
        (Default: `False`)
    over_clustering
        This argument can be provided in several ways:
        1) an input plain file with the over-clustering result of one cell per line.
        2) a string key specifying an existing metadata column in the AnnData (pre-created by the user).
        3) a python list, tuple, numpy array, pandas series or index representing the over-clustering result of the input cells.
        4) if none of the above is provided, will use a heuristic over-clustering approach according to the size of input data.
        Ignored if `majority_voting` is set to `False`.
    use_GPU
        Whether to use GPU for over clustering on the basis of `rapids-singlecell`. This argument is only relevant when `majority_voting = True`.
        (Default: `False`)
    min_prop
        For the dominant cell type within a subcluster, the minimum proportion of cells required to support naming of the subcluster by this cell type.
        Ignored if `majority_voting` is set to `False`.
        Subcluster that fails to pass this proportion threshold will be assigned `'Heterogeneous'`.
        (Default: 0)

    Returns
    ----------
    :class:`~celltypist.classifier.AnnotationResult`
        An :class:`~celltypist.classifier.AnnotationResult` object. Four important attributes within this class are:
        1) :attr:`~celltypist.classifier.AnnotationResult.predicted_labels`, predicted labels from celltypist.
        2) :attr:`~celltypist.classifier.AnnotationResult.decision_matrix`, decision matrix from celltypist.
        3) :attr:`~celltypist.classifier.AnnotationResult.probability_matrix`, probability matrix from celltypist.
        4) :attr:`~celltypist.classifier.AnnotationResult.adata`, AnnData representation of the input data.
    """
    #load model
    lr_classifier = model if isinstance(model, Model) else Model.load(model)
    #construct Classifier class
    clf = classifier.Classifier(filename = filename, model = lr_classifier, transpose = transpose_input, gene_file = gene_file, cell_file = cell_file)
    #predict
    predictions = clf.celltype(mode = mode, p_thres = p_thres)
    if not majority_voting:
        return predictions
    if predictions.cell_count <= 50:
        logger.warn(f"⚠️ Warning: the input number of cells ({predictions.cell_count}) is too few to conduct proper over-clustering; no majority voting is performed")
        return predictions
    #over clustering
    if over_clustering is None:
        over_clustering = clf.over_cluster(use_GPU = use_GPU)
        predictions.adata = clf.adata
    elif isinstance(over_clustering, str):
        if over_clustering in clf.adata.obs:
            over_clustering = clf.adata.obs[over_clustering]
        else:
            logger.info(f"👀 Did not identify '{over_clustering}' as a cell metadata column, assume it to be a plain text file")
            try:
                with open(over_clustering, 'rt') as f:
                    over_clustering = [x.strip() for x in f.readlines()]
            except Exception as e:
                raise Exception(
                        f"🛑 {e}")
    if len(over_clustering) != clf.adata.n_obs:
        raise ValueError(
                f"🛑 Length of `over_clustering` ({len(over_clustering)}) does not match the number of input cells ({clf.adata.n_obs})")
    #majority voting
    return classifier.Classifier.majority_vote(predictions, over_clustering, min_prop = min_prop)
