.. celltypist documentation master file, created by
   sphinx-quickstart on Thu Mar  4 16:45:59 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CellTypist's documentation!
======================================

.. include:: ../../README.md
   :parser: myst_parser.sphinx_
   :start-line: 2
   :end-line: 11

.. include:: ../../README.md
   :parser: myst_parser.sphinx_
   :start-line: 22
   :end-line: 32

.. include:: ../../CHANGELOG.md
   :parser: myst_parser.sphinx_

.. toctree::
   :maxdepth: 2
   :caption: Tutorials:
   :hidden:

   notebook/celltypist_tutorial
   notebook/celltypist_tutorial_ml
   notebook/celltypist_tutorial_cv
   notebook/celltypist_tutorial_harmonisation
   notebook/celltypist_tutorial_integration

.. toctree::
   :maxdepth: 2
   :caption: API (classification):
   :hidden:

   celltypist.train
   celltypist.annotate
   celltypist.dotplot
   celltypist.models.download_models
   celltypist.samples.downsample_adata

.. toctree::
   :maxdepth: 2
   :caption: API (harmonisation):
   :hidden:

   celltypist.harmonize
   celltypist.treeplot
   celltypist.sankeyplot

.. toctree::
   :maxdepth: 2
   :caption: API (integration):
   :hidden:

   celltypist.integrate

.. toctree::
   :maxdepth: 2
   :caption: Package organization:
   :hidden:

   celltypist.classifier.AnnotationResult
   celltypist.classifier.Classifier
   celltypist.models.Model
   celltypist.contro.align.DistanceAlignment
   celltypist.contro.distance.Distance
   celltypist.contro.pct.PredictiveClusteringTree
