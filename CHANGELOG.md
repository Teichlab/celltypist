# Changelog
*********************************
## CellTypist 1.7.0 (June 22, 2025)
- Add GPU option for the command line [#141](https://github.com/Teichlab/celltypist/issues/141)
- Add model subset function [#147](https://github.com/Teichlab/celltypist/issues/147)
- Avoid checking built-in models when model path is provided in command line [#151](https://github.com/Teichlab/celltypist/issues/151)
- Adjust leiden params to parallel scanpy
## CellTypist 1.6.3 (June 6, 2024)
- GPU support based on rapids-singlecell during the over-clustering step [#110](https://github.com/Teichlab/celltypist/issues/110)
- Load custom model directly [#60](https://github.com/Teichlab/celltypist/issues/60)
## CellTypist 1.6.2 (November 1, 2023)
- Reformat error and warning messages during prediction [#89](https://github.com/Teichlab/celltypist/issues/89)
- Fix cuml import during model training
- Fix categorical output from new versions of pandas during majority voting [#96](https://github.com/Teichlab/celltypist/issues/96)
## CellTypist 1.6.1 (September 25, 2023)
- Add cuML-based GPU support for model training [#80](https://github.com/Teichlab/celltypist/pull/80)
- Detect data input format only for the first 1000 cells
- Fix KeyError in [to_adata](https://celltypist.readthedocs.io/en/latest/celltypist.classifier.AnnotationResult.html#celltypist.classifier.AnnotationResult.to_adata) when setting `mode = 'prob match'` in [celltypist.annotate](https://celltypist.readthedocs.io/en/latest/celltypist.annotate.html)
- Add gene symbol to ID conversion file and expand README [#87](https://github.com/Teichlab/celltypist/issues/87)
## CellTypist 1.6.0 (August 5, 2023)
- Separate [CellHint](https://github.com/Teichlab/cellhint) from CellTypist
- Keep strict one-to-one orthologs during species model conversion
- Auto-fetch feature transfer files for model conversion
- Adapt code to sklearn>=1.3.0 [#75](https://github.com/Teichlab/celltypist/issues/75)
## CellTypist 1.5.3 (July 10, 2023)
- Fix DistanceMetric import for sklearn >= 1.3.0 [#73](https://github.com/Teichlab/celltypist/issues/73)
- Detect input format for [celltypist.dotplot](https://celltypist.readthedocs.io/en/latest/celltypist.dotplot.html)
- Require leidenalg >= 0.9.0
## CellTypist 1.5.2 (June 8, 2023)
- Patch log1p serialization issue and model download error
## CellTypist 1.5.1 (May 26, 2023)
- Fix error in not importing spmatrix in pct
## CellTypist 1.5.0 (May 16, 2023)
A more mature release before the formal CellTypist v2
## CellTypist 1.4.0 (May 2, 2023)
- Add modules of cell type harmonisation and integration (a release before the formal CellTypist v2)
## CellTypist 1.3.0 (December 14, 2022)
- Transfer numpy matrix to array during prediction for compatibility with sklearn >= 1.0 [#50](https://github.com/Teichlab/celltypist/issues/50)
## CellTypist 1.2.0 (August 22, 2022)
- Report all model url request errors (including timeout error) [#28](https://github.com/Teichlab/celltypist/issues/28)
- Add `with_mean` parameter in [celltypist.train](https://celltypist.readthedocs.io/en/latest/celltypist.train.html) to optimize RAM usage by not subtracting the mean [#29](https://github.com/Teichlab/celltypist/issues/29)
- Raise warning instead of error during prediction when the input expression matrix contains only a subset of genes
- Reformat [docs](https://celltypist.readthedocs.io/en/latest/?badge=latest#) with independent pages for main functions
- Support dot plot when the prediction comes from the 'prob match' mode [#33](https://github.com/Teichlab/celltypist/issues/33)
- Disable hard-coded `max_iter` in [celltypist.train](https://celltypist.readthedocs.io/en/latest/celltypist.train.html) with default values varying according to dataset sizes
## CellTypist 1.1.0 (June 09, 2022)
- Add citation (PMID: 35549406)
- Add Docker/Singularity container
- Fix the error of [celltypist.dotplot](https://celltypist.readthedocs.io/en/latest/celltypist.dotplot.html) when reference cluster is type of int [#22](https://github.com/Teichlab/celltypist/issues/22)
- Relax the criteria of detecting the plain file format
- Expand the docstring of [convert](https://celltypist.readthedocs.io/en/latest/celltypist.models.Model.html#celltypist.models.Model.convert) to explicitly support gene name transfer [#26](https://github.com/Teichlab/celltypist/issues/26)
- Automatically transpose (with warning) the input gene-by-cell matrix instead of raising an error
## CellTypist 1.0.0 (May 12, 2022)
CellTypist 1.0.0
## CellTypist 0.2.1 (Apr 24, 2022)
Final beta version
## CellTypist 0.2.0 (Feb 12, 2022)
Closer-to-mature version
## CellTypist 0.1.9 (Dec 14, 2021)
Close-to-mature version
