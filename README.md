<p align="left"><img src="https://github.com/Teichlab/celltypist/blob/main/docs/source/_static/img/logo_celltypist.png" width="250" height="85"></p>

[![Python Versions](https://img.shields.io/badge/python-3.6+-brightgreen.svg)](https://pypi.org/project/celltypist) [![Documentation Status](https://readthedocs.org/projects/celltypist/badge/?version=latest)](https://celltypist.readthedocs.io/en/latest/?badge=latest)

CellTypist is an automated cell type annotation tool for scRNA-seq datasets on the basis of logistic regression classifiers optimised by the stochastic gradient descent algorithm. CellTypist allows for cell prediction using either built-in (with a current focus on immune sub-populations) or custom models, in order to assist in the accurate classification of different cell types and subtypes.

# CellTypist website
Information of CellTypist can be also found in our CellTypist portal. [![Website www.celltypist.org](https://img.shields.io/website-up-down-brightgreen-red/http/shields.io.svg)](https://www.celltypist.org/)

# Interactive tutorials
[Using CellTypist for cell type classification ![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Teichlab/celltypist/blob/main/docs/notebook/celltypist_tutorial.ipynb)  
[Using CellTypist for multi-label classification ![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Teichlab/celltypist/blob/main/docs/notebook/celltypist_tutorial_ml.ipynb)  
[Best practice in large-scale cross-dataset label transfer using CellTypist ![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Teichlab/celltypist/blob/main/docs/notebook/celltypist_tutorial_cv.ipynb)

# Install CellTypist
### Using pip [![PyPI](https://img.shields.io/pypi/v/celltypist.svg?color=brightgreen&style=flat)](https://pypi.org/project/celltypist)
```console
pip install celltypist
```

### Using conda [![install with bioconda](https://img.shields.io/conda/vn/bioconda/celltypist.svg?color=brightgreen&style=flat)](https://anaconda.org/bioconda/celltypist)
```console
conda install -c bioconda -c conda-forge celltypist
```

# Usage (classification)

<details>
<summary><strong>1. Use in the Python environment</strong></summary>

+ <details>
  <summary><strong>1.1. Import the module</strong></summary>

  ```python
  import celltypist
  from celltypist import models
  ```
  </details>

+ <details>
  <summary><strong>1.2. Download available models</strong></summary>

  The models serve as the basis for cell type predictions. Information of available models can be also found [here](https://www.celltypist.org/models).
  ```python
  #Show all available models that can be downloaded and used.
  models.models_description()
  #Download a specific model, for example, `Immune_All_Low.pkl`.
  models.download_models(model = 'Immune_All_Low.pkl')
  #Download a list of models, for example, `Immune_All_Low.pkl` and `Immune_All_High.pkl`.
  models.download_models(model = ['Immune_All_Low.pkl', 'Immune_All_High.pkl'])
  #Update the models by re-downloading the latest versions if you think they may be outdated.
  models.download_models(model = ['Immune_All_Low.pkl', 'Immune_All_High.pkl'], force_update = True)
  #Show the local directory storing these models.
  models.models_path
  ```
  A simple way is to download all available models. Since each model is on average 1 megabyte (MB), we encourage the users to download all of them.
  ```python
  #Download all the available models.
  models.download_models()
  #Update all models by re-downloading the latest versions if you think they may be outdated.
  models.download_models(force_update = True)
  ```
  By default, a folder `.celltypist/` will be created in the user's home directory to store model files. A different path/folder can be specified by exporting the environment variable `CELLTYPIST_FOLDER` in your configuration file (e.g. in `~/.bash_profile`).
  ```bash
  #In the shell configuration file.
  export CELLTYPIST_FOLDER='/path/to/model/folder/'
  ```
  </details>

+ <details>
  <summary><strong>1.3. Overview of the models</strong></summary>

  All models are serialised in a binary format by [pickle](https://docs.python.org/3/library/pickle.html).
  ```python
  #Get an overview of the models that are downloaded in `1.2.`.
  #By default (`on_the_fly = False`), all possible models (even those that are not downloaded) are shown.
  models.models_description(on_the_fly = True)
  ```
  </details>

+ <details>
  <summary><strong>1.4. Inspect the model of interest</strong></summary>

  To take a look at a given model, load the model as an instance of the [Model](https://celltypist.readthedocs.io/en/latest/celltypist.models.Model.html) class as defined in CellTypist.
  ```python
  #Select the model from the above list. If the `model` argument is not provided, will default to `Immune_All_Low.pkl`.
  model = models.Model.load(model = 'Immune_All_Low.pkl')
  #The model summary information.
  model
  #Examine cell types contained in the model.
  model.cell_types
  #Examine genes/features contained in the model.
  model.features
  ```
  </details>

+ <details>
  <summary><strong>1.5. Celltyping based on the input of count table</strong></summary>

  CellTypist accepts the input data as a count table (cell-by-gene or gene-by-cell) in the format of `.txt`, `.csv`, `.tsv`, `.tab`, `.mtx` or `.mtx.gz`. A raw count matrix (reads or UMIs) is required. Non-expressed genes (if you are sure of their expression absence in your data) are suggested to be included in the input table as well, as they point to the negative transcriptomic signatures when compared with the model used.
  ```python
  #Get a demo test data. This is a UMI count csv file with cells as rows and gene symbols as columns.
  input_file = celltypist.samples.get_sample_csv()
  ```
  Assign the cell type labels from the model to the input test cells using the [celltypist.annotate](https://celltypist.readthedocs.io/en/latest/celltypist.annotate.html) function.
  ```python
  #Predict the identity of each input cell.
  predictions = celltypist.annotate(input_file, model = 'Immune_All_Low.pkl')
  #Alternatively, the model argument can be a previously loaded `Model` as in 1.4.
  predictions = celltypist.annotate(input_file, model = model)
  ```
  If your input file is in a gene-by-cell format (genes as rows and cells as columns), pass in the `transpose_input = True` argument. In addition, if the input is provided in the `.mtx` format, you will also need to specify the `gene_file` and `cell_file` arguments as the files containing names of genes and cells, respectively.
  ```python
  #In case your input file is a gene-by-cell table.
  predictions = celltypist.annotate(input_file, model = 'Immune_All_Low.pkl', transpose_input = True)
  #In case your input file is a gene-by-cell mtx file.
  predictions = celltypist.annotate(input_file, model = 'Immune_All_Low.pkl', transpose_input = True, gene_file = '/path/to/gene/file.txt', cell_file = '/path/to/cell/file.txt')
  ```
  Again, if the `model` argument is not specified, CellTypist will by default use the `Immune_All_Low.pkl` model.  
  
  The `annotate` function will return an instance of the [AnnotationResult](https://celltypist.readthedocs.io/en/latest/celltypist.classifier.AnnotationResult.html) class as defined in CellTypist.
  ```python
  #Summary information for the prediction result.
  predictions
  #Examine the predicted cell type labels.
  predictions.predicted_labels
  #Examine the matrix representing the decision score of each cell belonging to a given cell type.
  predictions.decision_matrix
  #Examine the matrix representing the probability each cell belongs to a given cell type (transformed from decision matrix by the sigmoid function).
  predictions.probability_matrix
  ```
  By default, with the `annotate` function, each query cell is predicted into the cell type with the largest score/probability among all possible cell types (`mode = 'best match'`). This mode is straightforward and can be used to differentiate between highly homogeneous cell types.  
  
  However, in some scenarios where a query cell cannot be assigned to any cell type in the reference model (i.e., a novel cell type) or can be assigned to multiple cell types (i.e., multi-label classification), a mode of probability match can be turned on (`mode = 'prob match'`) with a probability cutoff (default to 0.5, `p_thres = 0.5`) to decide the cell types (none, 1, or multiple) assigned for a given cell.
  ```python
  #Query cell will get the label of 'Unassigned' if it fails to pass the probability cutoff in each cell type.
  #Query cell will get multiple label outputs (concatenated by '|') if more than one cell type passes the probability cutoff.
  predictions = celltypist.annotate(input_file, model = 'Immune_All_Low.pkl', mode = 'prob match', p_thres = 0.5)
  ```
  The three tables in the `AnnotationResult` (`.predicted_labels`, `.decision_matrix` and `.probability_matrix`) can be written out to local files (tables) by the function [to_table](https://celltypist.readthedocs.io/en/latest/celltypist.classifier.AnnotationResult.html#celltypist.classifier.AnnotationResult.to_table), specifying the target `folder` for storage and the `prefix` common to each table.
  ```python
  #Export the three results to csv tables.
  predictions.to_table(folder = '/path/to/a/folder', prefix = '')
  #Alternatively, export the three results to a single Excel table (.xlsx).
  predictions.to_table(folder = '/path/to/a/folder', prefix = '', xlsx = True)
  ```
  The resulting `AnnotationResult` can be also transformed to an [AnnData](https://anndata.readthedocs.io/en/latest/) which stores the expression matrix in the log1p normalised format (to 10,000 counts per cell) by the function [to_adata](https://celltypist.readthedocs.io/en/latest/celltypist.classifier.AnnotationResult.html#celltypist.classifier.AnnotationResult.to_adata). The predicted cell type labels can be inserted to this `AnnData` as well by specifying `insert_labels = True` (which is the default behavior of `to_adata`).  
  
  Confidence scores of query cells can be inserted by specifying `insert_conf = True` (which is also the default behavior of `to_adata`). The scores correspond to the probabilities of cell predictions based on either `predictions.predicted_labels.predicted_labels` or `predictions.predicted_labels.majority_voting` (see `1.7.`), which can be specified by `insert_conf_by` (default to the former, `predicted_labels`).
  ```python
  #Get an `AnnData` with predicted labels and confidence scores embedded into the observation metadata columns.
  adata = predictions.to_adata(insert_labels = True, insert_conf = True)
  #Inspect these columns (`predicted_labels` and `conf_score`).
  adata.obs
  ```
  In addition, you can insert the decision matrix into the `AnnData` by passing in `insert_decision = True`, which represents the decision scores of each cell type distributed across the input cells. Alternatively, setting `insert_prob = True` will insert the probability matrix into the `AnnData`. The latter is the recommended way as probabilities are more interpretable (though sometimes not all query datasets converge to a meaningful range of probability values).  
    
  After the insertion, multiple columns will show up in the cell metadata of `AnnData`, with each column's name as a cell type name. Of note, all these columns (including the `predicted_labels` and `conf_score`) can be prefixed with a specific string by setting `prefix` in `to_adata`.
  ```python
  #Get an `AnnData` with predicted labels, confidence scores, and decision matrix.
  adata = predictions.to_adata(insert_labels = True, insert_conf = True, insert_decision = True)
  #Get an `AnnData` with predicted labels, confidence scores, and probability matrix (recommended).
  adata = predictions.to_adata(insert_labels = True, insert_conf = True, insert_prob = True)
  ```
  You can now manipulate this object with any functions or modules applicable to `AnnData`. Actually, CellTypist provides a quick function [to_plots](https://celltypist.readthedocs.io/en/latest/celltypist.classifier.AnnotationResult.html#celltypist.classifier.AnnotationResult.to_plots) to visualise your `AnnotationResult` and store the figures without the need of explicitly transforming it into an `AnnData`.
  ```python
  #Visualise the predicted cell types overlaid onto the UMAP.
  predictions.to_plots(folder = '/path/to/a/folder', prefix = '')
  ```
  A different prefix for the output figures can be specified with the `prefix` tag, and UMAP coordinates will be generated for the input dataset using a canonical [Scanpy](https://scanpy.readthedocs.io/en/stable/) pipeline. The labels in the figure may be crowded if too many cell types are predicted (can be alleviated by a majority voting process, see `1.7.`).  
    
  If you also would like to inspect the decision score and probability distributions for each cell type involved in the model, pass in the `plot_probability = True` argument. This may take a bit longer time as one figure will be generated for each of the cell types from the model.
  ```python
  #Visualise the decision scores and probabilities of each cell type overlaid onto the UMAP as well.
  predictions.to_plots(folder = '/path/to/a/folder', prefix = '', plot_probability = True)
  ```
  Multiple figures will be generated, including the predicted cell type labels overlaid onto the UMAP space, plus the decision score and probability distributions of each cell type on the UMAP.
  </details>

+ <details>
  <summary><strong>1.6. Celltyping based on AnnData</strong></summary>

  CellTypist also accepts the input data as an [AnnData](https://anndata.readthedocs.io/en/latest/) generated from for example [Scanpy](https://scanpy.readthedocs.io/en/stable/).  
    
  Since the expression of each gene will be centred and scaled by matching with the mean and standard deviation of that gene in the provided model, CellTypist requires a logarithmised and normalised expression matrix stored in the `AnnData` (log1p normalised expression to 10,000 counts per cell). CellTypist will try the `.X` attribute first, and if it does not suffice, try the `.raw.X` attribute. If none of them fit into the desired data type or the expression matrix is not properly normalised, an error will be raised.  
    
  Within the `AnnData`, please provide all genes to ensure maximal overlap with genes in the model. If you normalise and logarithmise the gene expression matrix using all genes while later only keep a subset of genes in the `AnnData`, the prediction result may not be optimal.
  ```python
  #Provide the input as an `AnnData`.
  predictions = celltypist.annotate('/path/to/input.h5ad', model = 'Immune_All_Low.pkl')
  #Alternatively, the input can be specified as an `AnnData` already loaded in memory.
  predictions = celltypist.annotate(a_loaded_adata, model = 'Immune_All_Low.pkl')
  ```
  All the parameters and downstream operations are the same as in `1.5.`, except that 1) the transformed `AnnData` from `to_adata` stores all the expression matrix and other information as is in the original object. 2) when generating the visualisation figures, existing UMAP coordinates will be used. If no UMAP coordinates are found, CellTypist will fall back on the neighborhood graph to yield new 2D UMAP projections. If none is available, a canonical Scanpy pipeline will be performed to generate the UMAP coordinates as in `1.5.`.  
    
  Of note, when the input is an `AnnData`, compared to the visualisations in `1.5.`, a more useful way for visualising the prediction result is to use the function `celltypist.dotplot`, which quantitatively compares the CellTypist prediction result with the cell types (or clusters) pre-defined and stashed in the `AnnData` by the user. Specifically, a dot plot will be generated, demonstrating the match between CellTypist predictions and manual annotations (or clusters). For each cell type or cluster (each column within the dot plot), this plot shows how it can be 'decomposed' into different cell types predicted by CellTypist.
  ```python
  #Examine the correspondence between CellTypist predictions (`use_as_prediction`) and manual annotations (`use_as_reference`).
  #Here, `predicted_labels` from `predictions.predicted_labels` is used as the prediction result from CellTypist.
  #`use_as_prediction` can be also set as `majority_voting` (see `1.7.`).
  celltypist.dotplot(predictions, use_as_reference = 'column_key_of_manual_annotation', use_as_prediction = 'predicted_labels')
  ```
  Check [celltypist.dotplot](https://celltypist.readthedocs.io/en/latest/celltypist.dotplot.html) for other parameters controlling visualisation details of this plot.
  </details>

+ <details>
  <summary><strong>1.7. Use a majority voting classifier combined with celltyping</strong></summary>

  By default, CellTypist will only do the prediction jobs to infer the identities of input cells, which renders the prediction of each cell independent. To combine the cell type predictions with the cell-cell transcriptomic relationships, CellTypist offers a majority voting approach based on the idea that similar cell subtypes are more likely to form a (sub)cluster regardless of their individual prediction outcomes.
  To turn on the majority voting classifier in addition to the CellTypist predictions, pass in `majority_voting = True` to the `annotate` function.
  ```python
  #Turn on the majority voting classifier as well.
  predictions = celltypist.annotate(input_file, model = 'Immune_All_Low.pkl', majority_voting = True)
  ```
  During the majority voting, to define cell-cell relations, CellTypist will use a heuristic over-clustering approach according to the size of the input data with the aid of a Leiden clustering pipeline. Users can also provide their own over-clustering result to the `over_clustering` argument. This argument can be specified in several ways:
   1) an input plain file with the over-clustering result of one cell per line.
   2) a string key specifying an existing cell metadata column in the `AnnData` (pre-created by the user).
   3) a list-like object (such as a numpy 1D array) indicating the over-clustering result of all cells.
   4) if none of the above is provided, will use a heuristic over-clustering approach, noted above.
  ```python
  #Add your own over-clustering result.
  predictions = celltypist.annotate(input_file, model = 'Immune_All_Low.pkl', majority_voting = True, over_clustering = '/path/to/over_clustering/file')
  ```
  There is also a `min_prop` parameter (defaults to 0) which controls the minimum proportion of cells from the dominant cell type required to name a given subcluster by this cell type. Subcluster that fails to pass this proportion threshold will be assigned `Heterogeneous`.  
    
  Similarly, an instance of the `AnnotationResult` class will be returned.
  ```python
  #Examine the predicted cell type labels.
  predictions.predicted_labels
  #Examine specifically the majority-voting results.
  predictions.predicted_labels.majority_voting
  #Examine the matrix representing the decision score of each cell belonging to a given cell type.
  predictions.decision_matrix
  #Examine the matrix representing the probability each cell belongs to a given cell type (transformed from decision matrix by the sigmoid function).
  predictions.probability_matrix
  ```
  Compared to the results without majority-voting functionality as in `1.5.` and `1.6.`, the `.predicted_labels` attribute now has two extra columns (`over_clustering` and `majority_voting`) in addition to the column `predicted_labels`.  
    
  Other parameters and downstream operations are the same as in `1.5.` and `1.6.`. Note that due to the majority-voting results added, the exported tables (by `to_table`), the transformed `AnnData` (by `to_adata`), and the visualisation figures (by `to_plots`) will all have additional outputs or information indicating the majority-voting outcomes. For example, when using the function `celltypist.dotplot`, you can set `use_as_prediction = 'majority_voting'` to visualise the match between majority-voting results with manual annotations. The other example is that when using `to_adata`, you can specify `insert_conf_by = 'majority_voting'` to have the confidence scores corresponding to the majority-voting result instead of raw predictions (`insert_conf_by = 'predicted_labels'` which is the default).
  ```python
  #Examine the correspondence between CellTypist predictions (`use_as_prediction`) and manual annotations (`use_as_reference`).
  celltypist.dotplot(predictions, use_as_reference = 'column_key_of_manual_annotation', use_as_prediction = 'majority_voting')
  ```
  </details>
</details>

<details>
<summary><strong>2. Use as the command line</strong></summary>

+ <details>
  <summary><strong>2.1. Check the command line options</strong></summary>

  ```bash
  celltypist --help
  ```
  </details>

+ <details>
  <summary><strong>2.2. Download all available models</strong></summary>

  ```bash
  celltypist --update-models
  ```
  This will download the latest models from the remote server.
  </details>

+ <details>
  <summary><strong>2.3. Overview of the models</strong></summary>

  ```bash
  celltypist --show-models
  ```
  </details>

+ <details>
  <summary><strong>2.4. Celltyping based on the input of count table</strong></summary>

  See `1.5.` for the format of the desired count matrix.
  ```bash
  celltypist --indata /path/to/input/file --model Immune_All_Low.pkl --outdir /path/to/outdir
  ```
  You can add a different model to be used in the `--model` option. If the `--model` is not provided, CellTypist will by default use the `Immune_All_Low.pkl` model. The output directory will be set to the current working directory if `--outdir` is not specified.  
    
  If your input file is in a gene-by-cell format (genes as rows and cells as columns), add the `--transpose-input` option.
  ```bash
  celltypist --indata /path/to/input/file --model Immune_All_Low.pkl --outdir /path/to/outdir --transpose-input
  ```
  If the input is provided in the `.mtx` format, you will also need to specify the `--gene-file` and `--cell-file` options as the files containing names of genes and cells, respectively.  
    
  The default mode (`--mode best_match`) for prediction is to choose the cell type with the largest score/probability as the final prediction; setting `--mode prob_match` combined with a probability threshold (default to 0.5, `--p-thres 0.5`) will enable a multi-label classification, which assigns 0 (i.e., unassigned), 1, or >=2 cell type labels to each query cell.  
    
  Other options that control the output files of CellTypist include `--prefix` which adds a custom prefix and `--xlsx` which merges the output files into one xlsx table. Check `celltypist --help` for more details.
  </details>

+ <details>
  <summary><strong>2.5. Celltyping based on AnnData</strong></summary>

  See `1.6.` for the requirement of the expression matrix in the AnnData object (`.h5ad`).
  ```bash
  celltypist --indata /path/to/input/adata --model Immune_All_Low.pkl --outdir /path/to/outdir
  ```
  Other command line options are the same as in `2.4.`.
  </details>

+ <details>
  <summary><strong>2.6. Use a majority voting classifier combined with celltyping</strong></summary>

  See `1.7.` for how the majority voting classifier works.
  ```bash
  celltypist --indata /path/to/input/file --model Immune_All_Low.pkl --outdir /path/to/outdir --majority-voting
  ```
  During the majority voting, to define cell-cell relations, CellTypist will use a heuristic over-clustering approach according to the size of the input data with the aid of a Leiden clustering pipeline. Users can also provide their own over-clustering result to the `--over-clustering` option. This option can be specified in several ways:
     1) an input plain file with the over-clustering result of one cell per line.
     2) a string key specifying an existing cell metadata column in the `AnnData` (pre-created by the user).
     3) if none of the above is provided, will use a heuristic over-clustering approach, noted above.
  ```bash
  celltypist --indata /path/to/input/file --model Immune_All_Low.pkl --outdir /path/to/outdir --majority-voting --over-clustering /path/to/over_clustering/file
  ```
  There is also a `--min-prop` option (defaults to 0) which controls the minimum proportion of cells from the dominant cell type required to name a given subcluster by this cell type. Subcluster that fails to pass this proportion threshold will be assigned `Heterogeneous`.  
    
  Other command line options are the same as in `2.4.`.
  </details>

+ <details>
  <summary><strong>2.7. Generate visualisation figures for the results</strong></summary>

  In addition to the tables output by CellTypist, you have the option to generate multiple figures to get an overview of your prediction results. See `1.5.`, `1.6.` and `1.7.` for what these figures represent.
  ```bash
  #Plot the results after the celltyping process.
  celltypist --indata /path/to/input/file --model Immune_All_Low.pkl --outdir /path/to/outdir --plot-results
  #Plot the results after the celltyping and majority-voting processes.
  celltypist --indata /path/to/input/file --model Immune_All_Low.pkl --outdir /path/to/outdir --majority-voting --plot-results
  ```
  </details>
</details>

<details>
<summary><strong>3. Use in the R environment</strong></summary>

Currently, there is no plan for R compatibility. Try to convert R objects into AnnData for use in CellTypist.
</details>

<details>
<summary><strong>4. Use as Docker/Singularity container</strong></summary>

  ### Docker

  A docker image is available from the Quay.io Container Registry as [`quay.io/teichlab/celltypist:latest`](https://quay.io/repository/teichlab/celltypist?tab=tags).
  
  **Simple usage:**
  ```bash
  docker run --rm -it \
    -v /path/to/data:/data \
    quay.io/teichlab/celltypist:latest \
    celltypist --indata /data/file --model Immune_All_Low.pkl --outdir /data/output
  ```
  **Usage with custom models:**
  ```bash
  docker run --rm -it \
    -v /path/to/data:/data \
    -v /path/to/models:/opt/celltypist/data/models \
    quay.io/teichlab/celltypist:latest \
    celltypist --indata /data/file --model My_Custom_Model.pkl --outdir /data/output
  ```
  
  ### Singularity
  
  Use the `singularity pull` command to download the container from the given container registry:
  ```bash
  singularity pull celltypist-latest.sif docker://quay.io/teichlab/celltypist:latest
  ```
  Then run the downloaded image as a container.
  
  **Simple usage:**
  ```bash
  singularity run \
    -B /path/to/data:/data \
    celltypist-latest.sif \
    celltypist --indata /data/file --model Immune_All_Low.pkl --outdir /data/output
  ```
  **Usage with custom models:**
  ```bash
  singularity run \
    -B /path/to/data:/data \
    -B /path/to/models:/opt/celltypist/data/models \
    celltypist-latest.sif \
    celltypist --indata /data/file --model My_Custom_Model.pkl --outdir /data/output
  ```
  
</details>

<details>
<summary><strong>Supplemental guidance</strong></summary>

+ <details>
  <summary><strong>Generate a custom model</strong></summary>
  
  As well as the models provided by CellTypist (see `1.2.`), you can generate your own model from which the cell type labels can be transferred to another scRNA-seq dataset. This will be most useful when a large and comprehensive reference atlas is trained for future use, or when the similarity between two scRNA-seq datasets is under examination.  
    
  ### Inputs for data training
  The inputs for CellTypist training comprise the gene expression data, the cell annotation details (i.e., cell type labels), and in some scenarios the genes used. To facilitate the training process, the `train` function (see below) has been designed to accommodate different kinds of input formats:
     1) The gene expression data can be provided as a path to the expression table (such as `.csv` and `.mtx`), or a path to the `AnnData` (`.h5ad`), with the former containing raw counts (in order to reduce the file size) while the latter containing log1p normalised expression (to 10,000 counts per cell) stored in `.X` or `.raw.X`. In addition to specifying the paths, you can provide any array-like objects (e.g., `csr_matrix`) or `AnnData` which are already loaded in memory (both should be in the log1p format). A cell-by-gene format (cells as rows and genes as columns) is required.
     2) The cell type labels can be supplied as a path to the file containing cell type label per line corresponding to the cells in gene expression data. Any list-like objects (such as a `tuple` or `series`) are also acceptable. If the gene expression data is input as an `AnnData`, you can also provide a column name from its cell metadata (`.obs`) which represents information of cell type labels.
     3) The genes will be automatically extracted if the gene expression data is provided as a table file, an `AnnData` or a `DataFrame`. Otherwise, you need to specify a path to the file containing one gene per line corresponding to the genes in the gene expression data. Any list-like objects (such as a `tuple` or `series`) are also acceptable.
  
  ### One-pass data training
  Derive a new model by training the data using the [celltypist.train](https://celltypist.readthedocs.io/en/latest/celltypist.train.html) function:
  ```python
  #Training a CellTypist model.
  new_model = celltypist.train(expression_input, labels = label_input, genes = gene_input)
  ```
  If the input is a table file, an `AnnData` or a `DataFrame`, genes will be automatically extracted and the `genes` tag can thus be omitted from the above code. If your input is in a gene-by-cell format (genes as rows and cells as columns), remember to pass in the `transpose_input = True` argument.  
    
  Before the training is conducted, the gene expression format will be checked to make sure the input data is supplied as required. For example, the expression matrix should be in log1p normalised expression (to 10,000 counts per cell) if the input is an `AnnData`. This means when you subset the input with given genes (e.g., by highly variable genes), an error may be raised as CellTypist cannot judge the input as properly normalised with only a subset of genes. In such a case, pass in `check_expression = False` to skip the expression format check.
  ```python
  #Training a CellTypist model with only subset of genes (e.g., highly variable genes).
  #Restricting the input to a subset of genes can accelerate the training process.
  #Use `AnnData` here as an example.
  new_model = celltypist.train(some_adata[:, some_adata.var.highly_variable], labels = label_input, check_expression = False)
  ```
  By default, data is trained using a traditional logistic regression classifier. This classifier is well suited to datasets of small or intermediate sizes (as an empirical estimate, <= 100k cells), and usually leads to an unbiased probability range with less parameter tuning. Among the training parameters, three important ones are `solver` which (if not specified by the user) is selected based on the size of the input data by CellTypist, `C` which sets the inverse of L2 regularisation strength, and `max_iter` which controls the maximum number of iterations before reaching the minimum of the cost function. Other (hyper)parameters from [LogisticRegression](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html) are also applicable in the `train` function.  
    
  When the dimensions of the input data are large, training may take longer time even with CPU parallelisation (achieved by the `n_jobs` argument). To reduce the training time as well as to add some randomness to the classifier's solution, a stochastic gradient descent (SGD) logistic regression classifier can be enabled by `use_SGD = True`.
  ```python
  #Training a CellTypist model with SGD learning.
  new_model = celltypist.train(expression_input, labels = label_input, genes = gene_input, use_SGD = True)
  ```
  A logistic regression classifier with SGD learning reduces the training burden dramatically and has a comparable performance versus a traditional logistic regression classifier. A minor caveat is that more careful model parameter tuning may be needed if you want to utilise the probability values from the model for scoring cell types in the prediction step (the selection of the most likely cell type for each query cell is not influenced however). Among the training parameters, two important ones are `alpha` which sets the L2 regularisation strength and `max_iter` which controls the maximum number of iterations. Other (hyper)parameters from [SGDClassifier](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.SGDClassifier.html) are also applicable in the `train` function.  
    
  When the training data contains a huge number of cells (for example >500k cells) or more randomness in selecting cells for training is needed, you may consider using the mini-batch version of the SGD logistic regression classifier by specifying `use_SGD = True` and `mini_batch = True`. As a result, in each epoch (default to 10 epochs, `epochs = 10`), cells are binned into equal-sized (the size is default to 1000, `batch_size = 1000`) random batches, and are trained in a batch-by-batch manner (default to 100 batches, `batch_number = 100`).
  ```python
  #Get a CellTypist model with SGD mini-batch training.
  new_model = celltypist.train(expression_input, labels = label_input, genes = gene_input, use_SGD = True, mini_batch = True)
  ```
  By selecting part of cells for training (default to 1,000,000 cells with possible duplications, `epochs` x `batch_size` x `batch_number`), training time can be again reduced and the performance of the derived model is shown to persist as compared to the above two methods. Since some rare cell types may be undersampled during this procedure, you can pass in the `balance_cell_type = True` argument to sample rare cell types with a higher probability, ensuring close-to-even cell type distributions in mini-batches (subject to the maximum number of cells that can be provided by a given cell type).
    
  There are also some free texts that can be inserted (e.g., `date`) to describe the model. Check out the [celltypist.train](https://celltypist.readthedocs.io/en/latest/celltypist.train.html) for more information.  
    
  The resulting model is an instance of the `Model` class as in `1.4.`, and can be manipulated as with other CellTypist models.  
    
  Save this model locally:
  ```python
  #Write out the model.
  new_model.write('/path/to/local/folder/some_model_name.pkl')
  ```
  A suggested location for stashing the model is the `models.models_path` (see `1.2.`). Through this, all models (including the models provided by CellTypist) will be in the same folder, and can be accessed in the same manner as in `1.4.`.
  ```python
  #Write out the model in the `models.models_path` folder.
  new_model.write(f'{models.models_path}/some_model_name.pkl')
  ```
  To leverage this model, first load it by `models.Model.load`.
  ```python
  new_model = models.Model.load('/path/to/local/folder/some_model_name.pkl')
  ```
  This model can be used as with the built-in CellTypist models, for example, it can be specified as the `model` argument in `annotate`.
  ```python
  #Predict the identity of each input cell with the new model.
  predictions = celltypist.annotate(input_file, model = new_model)
  #Alternatively, just specify the model path (recommended as this ensures the model is intact every time it is loaded).
  predictions = celltypist.annotate(input_file, model = '/path/to/local/folder/some_model_name.pkl')
  #If the model is stored in `models.models_path`, only the model name is needed.
  predictions = celltypist.annotate(input_file, model = 'some_model_name.pkl')
  ```
  Downstream operations are the same as in `1.4.`, `1.5.`, `1.6.`, and `1.7.`.
  
  ### Two-pass data training incorporating feature selection
  Some scRNA-seq datasets may involve the noise mostly from genes not helpful or even detrimental to the characterisation of cell types. To mitigate this, `celltypist.train` has the option (`feature_selection = True`) to do a fast feature selection based on the feature importance (here, the absolute regression coefficients) using SGD learning. In short, top important genes (default: `top_genes = 300`) are selected from each cell type, and are further combined across cell types as the final feature set. The classifier is then re-run using the corresponding subset of the input data.
  ```python
  #Two-pass data training with traditional logistic regression after SGD-based feature selection.
  new_model = celltypist.train(expression_input, labels = label_input, genes = gene_input, feature_selection = True)
  #Two-pass data training with SGD learning after feature selection.
  new_model = celltypist.train(expression_input, labels = label_input, genes = gene_input, use_SGD = True, feature_selection = True)
  #Two-pass data training with SGD mini-batch training after feature selection.
  new_model = celltypist.train(expression_input, labels = label_input, genes = gene_input, use_SGD = True, mini_batch = True, feature_selection = True)
  ```
  If you prefer other feature selection approaches and obtain a set of genes which are designated as important features, you can subset your input data and train the CellTypist model accordingly. As noted in the previous section, remember to pass in the `check_expression = False` argument.
  ```python
  new_model = celltypist.train(expression_input_subset, labels = label_input, genes = gene_input, check_expression = False)
  ```
  The downstream workflow is the same as that from one-pass data training.

  ### General parameters relating to runtime and RAM usage
  `max_iter`: when `celltypist.train` does not converge for a long time, setting `max_iter` to a lower number can reduce runtime at a possible cost of a suboptimal model.  
    
  `with_mean`: when the training data is a sparse matrix, setting `with_mean = False` will preserve sparsity by skipping the step of subtraction by the mean during scaling, and thus lower the RAM usage at the cost of a suboptimal model.  
    
  `n_jobs`: Number of CPUs used. This argument is not applicable to mini-batch training.  
    
  `use_GPU`: GPU acceleration by using logistic regression from [cuml](https://docs.rapids.ai/api/cuml/stable). You need to install RAPIDS and cuml first. This argument is ignored if SGD learning is enabled.
  </details>

+ <details>
  <summary><strong>Cross-species model conversion</strong></summary>

  It is always recommended to predict a query dataset using the reference model from the same species. In cases where a cross-species label projection is needed, you can convert the model of interest to its "orthologous" form of another species. This is achieved by aligning orthologous genes between species.  
    
  Load a human immune model.
  ```python
  model = models.Model.load('Immune_All_Low.pkl')
  ```
  This model can be converted to a mouse equivalent through the [convert](https://celltypist.readthedocs.io/en/latest/celltypist.models.Model.html#celltypist.models.Model.convert) method. By default, a human-mouse conversion (or the opposite) will be conducted by automatically detecting the species of the model (e.g., human) and transforming it to the other species (e.g., mouse).
  ```python
  #Note `model` is modified in-place.
  model.convert()
  ```
  By default (`unique_only = True`), only 1:1 orthologs between the two species are kept and all other genes are discarded in the model. You can also keep those genes (including both 1:N and N:1 orthologs) by specifying `unique_only = False`. By doing so, you need to specify how these 1:N orthologs will be handled: for each gene, averaging the classifier weights (`collapse = 'average'`, which is the default when `unique_only = False`) or randomly choosing one gene's weight as the representative (`collapse = 'random'`) from all its orthologs.
  ```python
  #For illustration purpose. Convert the model by utilising 1:N orthologs and their average weights.
  #model.convert(unique_only = False, collapse = 'average')
  ```
  As mentioned above, the default mode is a human-to-mouse (or mouse-to-human) conversion using the built-in gene mapping [file](https://github.com/Teichlab/celltypist/blob/main/celltypist/data/samples/Ensembl105_Human2Mouse_Genes.csv) (Ensembl105 version). For conversion to other species, you can provide a different file (`map_file`), with one column being the species of the model and the other column being the species you want to convert to. Check out `models.Model.convert` for more information.  
    
  Lastly, write out the converted model locally.
  ```python
  model.write('/path/to/local/folder/some_model_name.pkl')
  ```
  This model can be used as with other CellTypist models.
  </details>

+ <details>
  <summary><strong>Model conversion from gene symbols to Ensembl IDs</strong></summary>

  CellTypist models are usually trained based on gene symbols. When genes of a query dataset are formatted as Ensembl IDs, you can convert gene symbols in the model to Ensembl ID for matching the query dataset. The [convert](https://celltypist.readthedocs.io/en/latest/celltypist.models.Model.html#celltypist.models.Model.convert) method will be utilised as in the above section.  
    
  Specifically, you need to provide a gene-symbol-to-Ensembl-ID file, such that gene symbols in the model will be converted to IDs (or vice versa). A built-in [file](https://github.com/Teichlab/celltypist/blob/main/celltypist/data/samples/GENCODEv44_Gene_id2name.csv) is provided in CellTypist (GENCODE v44). Parameters and details during model conversion can be found in the previous section `Cross-species model conversion`.  
    
  Load a human immune model.
  ```python
  model = models.Model.load('Immune_All_Low.pkl')
  ```
  Convert gene symbols to Ensembl IDs using the built-in file. You can also provide a path to your own ID mapping file.
  ```python
  #Note `model` is modified in-place.
  model.convert('GENCODEv44_Gene_id2name.csv')
  ```
  Lastly, write out the converted model locally.
  ```python
  model.write('/path/to/local/folder/some_model_name.pkl')
  ```
  This model can be used as with other CellTypist models.
  </details>
</details>

# Citation
Dominguez Conde et al., Cross-tissue immune cell analysis reveals tissue-specific features in humans. Science 376, eabl5197 (2022). [Link](https://doi.org/10.1126/science.abl5197)
