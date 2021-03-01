from . import classifier, models

def annotate(filename: str,
             model: str = "",
             transpose_input: bool = False,
             majority_voting: bool = False,
             over_clustering = None) -> classifier.AnnotationResult:
    """Run the prediction and (optional) majority voting to annotate the input dataset."""
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
