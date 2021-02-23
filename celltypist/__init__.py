from celltypist import classifier
from celltypist import models
from celltypist import defaults
from celltypist import samples
from celltypist import logger
from celltypist import helpers

__version__ = "0.1.9"

def annotate(filename: str,
             model: str = "", transpose_input = False):
             #,
             #chunk_size: int = defaults.chunk_size,
             #cpus: int = defaults.max_cpus,
             #quiet: bool = False) -> classifier.AnnotationResult:
    """Run celltying process to annotate the dataset."""
    sgd_classifier = models.load(model)

    clf = classifier.Classifier(
        filename=filename,
        model=sgd_classifier, transpose = transpose_input)
        # ,
        # cpus=cpus,
        # chunk_size=chunk_size,
        # quiet=quiet)

    return clf.celltype()
