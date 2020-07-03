from celltypist import classifier
from celltypist import models
from celltypist import defaults
from celltypist import samples


def annotate(filename: str,
             model: str = defaults.model,
             chunk_size: int = defaults.chunk_size,
             cpus: int = defaults.max_cpus) -> classifier.AnnotationResult:
    """Run celltying process to annotate the dataset."""
    sgd_classifier = models.load(model)
    annotator = classifier.Classifier(
        filename, model=sgd_classifier, cpus=min(cpus, defaults.max_cpus),
        chunk_size=chunk_size)
    return annotator.celltype()
