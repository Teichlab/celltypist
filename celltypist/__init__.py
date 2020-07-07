from celltypist import classifier
from celltypist import models
from celltypist import defaults
from celltypist import samples
from celltypist import logger


def annotate(filename: str,
             model: str = defaults.model,
             chunk_size: int = defaults.chunk_size,
             cpus: int = defaults.max_cpus,
             quiet: bool = False) -> classifier.AnnotationResult:
    """Run celltying process to annotate the dataset."""
    sgd_classifier = models.load(model)

    clf = classifier.Classifier(
        filename=filename,
        model=sgd_classifier,
        cpus=cpus,
        chunk_size=chunk_size,
        quiet=quiet)

    return clf.celltype()
