from . import classifier, models

def annotate(filename: str,
             model: str = "",
             transpose_input: bool = False,
             majority_voting: bool = False,
             over_clustering = None) -> classifier.AnnotationResult:
    """
    Run the prediction and (optional) majority voting to annotate the input dataset.

    Parameters
    ----------
    filename
        Path to the input count matrix (supported types are csv, txt, tsv and tab) or Scanpy object (h5ad).
        If it's the former, a cell-by-gene format is desirable (see `transpose_input` for more information).
        Genes should be gene symbols. Non-expressed genes are preferred to be provided as well.
    model
        Model used to predict the input cells. Default to using the `Immune_All_Low.pkl` model.
        To see all available models and their descriptions, use `celltypist.models.models_description()`.
    transpose_input
        Whether to transpose the input matrix. Set to `True` if `filename` is provided in a gene-by-cell format.
        (Default: `False`)
    majority_voting
        Whether to refine the predicted labels by running the majority voting classifier after over-clustering.
        (Default: `False`)
    over_clustering
        This argument can be provided in several ways:
        1) an input plain file with the over-clustering result of one cell per line.
        2) a string key specifying an existing metadata column in the AnnData (pre-created by the user).
        3) a python list, numpy array, or pandas series indicating the over-clustering result of all cells.
        4) if none of the above is provided, will use a heuristic over-clustering approach based on input data size.
        Ignored if `majority_voting` is set to `False`.

    Returns
    ----------
    A `~celltypist.classifier.AnnotationResult` object. Two important attributes within are:
        1) `.predicted_labels`: predicted labels from celltypist.
        2) `.probability_table`: probability matrix from celltypist.
    """
    #load model
    sgd_classifier = models.load(model)
    #construct Classifier class
    clf = classifier.Classifier(filename = filename, model = sgd_classifier, transpose = transpose_input)
    #predict
    predictions = clf.celltype()
    if majority_voting is False:
        return predictions
    #over clustering
    if over_clustering is None:
        over_clustering = clf.over_cluster()
    elif len(over_clustering) == 1:
        if over_clustering in clf.adata.obs.columns.tolist():
            over_clustering = clf.adata.obs[over_clustering]
        else:
            try:
                with open(over_clustering, 'rt') as f:
                    over_clustering = [x.strip() for x in f.readlines()]
            except Exception as e:
                raise Exception(f"🛑 {e}")
    if len(over_clustering) != clf.adata.shape[0]:
        raise ValueError(f"🛑 length of over_clustering ({len(over_clustering)}) does not match the number of input cells ({clf.adata.shape[0]})")
    #majority voting
    return classifier.Classifier.majority_vote(predictions, over_clustering)
