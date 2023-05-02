# Changelog
*********************************
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
