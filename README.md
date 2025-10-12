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
  A simple way is to download all available models. Since each model is on average several megabytes (MB), we encourage the users to download all of them.
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
     2) The cell type labels can be supplied as a path to the file containing cell type label per line corresponding to the cells in gene expression data. Any list-like objects (such as a `tuple` or `Series`) are also acceptable. If the gene expression data is input as an `AnnData`, you can also provide a column name from its cell metadata (`.obs`) which represents information of cell type labels.
     3) The genes will be automatically extracted if the gene expression data is provided as a table file, an `AnnData` or a `DataFrame`. Otherwise, you need to specify a path to the file containing one gene per line corresponding to the genes in the gene expression data. Any list-like objects (such as a `tuple` or `Series`) are also acceptable.
  
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

+ <details>
  <summary><strong>Subset a model</strong></summary>

  A CellTypist model can be restricted to a specified subset of cell types, allowing you to limit the search space and exclude cell types unlikely to be present in your query data. This is achieved by the [subset](https://celltypist.readthedocs.io/en/latest/celltypist.models.Model.html#celltypist.models.Model.subset) method.  
    
  Note that this strategy is suboptimal/not recommended for subsetting a model. A more accurate approach is to retrain the original reference data using only the desired subset of cell types by [celltypist.train](https://celltypist.readthedocs.io/en/latest/celltypist.train.html).  
    
  Load a human immune model.
  ```python
  model = models.Model.load('Immune_All_Low.pkl')
  ```
  Keep only T cells in the model.
  ```python
  #Note `model` is modified in-place.
  model.subset(keep_cell_types = [x for x in model.cell_types if 'T cells' in x])
  ```
  Or alternatively, you can choose to exclude T cells from the model.
  ```python
  model.subset(exclude_cell_types = [x for x in model.cell_types if 'T cells' in x])
  ```
  Lastly, write out the sub-model locally.
  ```python
  model.write('/path/to/local/folder/some_model_name.pkl')
  ```
  This sub-model can be used as with other CellTypist models.
  </details>
</details>

# Usage (hierarchical classification)

<details>
<summary><strong>1. Use in the Python environment</strong></summary>

+ <details>
  <summary><strong>1.1. Import the module</strong></summary>

  ```python
  import celltypist
  from celltypist import models
  from celltypist.models import Model, HierModel
  ```
  </details>

+ <details>
  <summary><strong>1.2. Download available models</strong></summary>

  The models serve as the basis for cell type predictions. Information of available models can be also found [here](https://www.celltypist.org/models).  
  
  There are two types of models: `flat` and `hierarchical`, corresponding to the flat models in CellTypist v1.0 and the hierarchical models in CellTypist v2.0. Flat models assume independence among cell type labels or annotations, whereas hierarchical models encode predefined relationships among cell types within a structured [Tree](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html) (for details on what a tree is and how to construct, modify, and visualise it, see `4.`).  

  This section focuses on the use of built-in hierarchical models (for advanced topics such as inspecting, modifying, and training a hierarchical model, see `5.`). In brief, a hierarchical model is an ensemble model composed of multiple flat models, operating in one of two modes: Local Classifier per Parent Node (`LCPN`), which trains a local classifier for each parent node (coarse cell type), and Local Classifier per Level (`LCL`), which trains a local classifier at each level (depth) of the tree.
  ```python
  #Show all available models that can be downloaded and used.
  models.models_description()
  #Download a specific model, for example, `Human_Tissue_Immune_LCPN.pkl`.
  models.download_models(model = 'Human_Tissue_Immune_LCPN.pkl')
  #Download a list of models, for example, `Human_Tissue_Immune_LCPN.pkl` and `Human_Tissue_Immune_LCL.pkl`.
  models.download_models(model = ['Human_Tissue_Immune_LCPN.pkl', 'Human_Tissue_Immune_LCL.pkl'])
  #Update the models by re-downloading the latest versions if you think they may be outdated.
  models.download_models(model = ['Human_Tissue_Immune_LCPN.pkl', 'Human_Tissue_Immune_LCL.pkl'], force_update = True)
  #Show the local directory storing these models.
  models.models_path
  ```
  A simple way is to download all available models. Since each model is on average several megabytes (MB), we encourage the users to download all of them.
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
  <summary><strong>1.3. Load the model of interest</strong></summary>

  All models are serialised in a binary format by [pickle](https://docs.python.org/3/library/pickle.html).
  ```python
  #Get an overview of the models that are downloaded in `1.2.`.
  #By default (`on_the_fly = False`), all possible models (even those that are not downloaded) are shown.
  models.models_description(on_the_fly = True)
  ```
  To take a look at a given hierarchical model, load the model as an instance of the [HierModel](https://celltypist.readthedocs.io/en/latest/celltypist.models.HierModel.html) class as defined in CellTypist.
  ```python
  #Select the model from the above list (type = 'hierarchical'). If the `model` argument is not provided, will default to `Human_Tissue_Immune_LCPN.pkl`.
  hier_model = HierModel.load(model = 'Human_Tissue_Immune_LCPN.pkl')
  #The model summary information.
  hier_model
  #Examine the mode ('LCPN' or 'LCL') of the model.
  hier_model.mode
  #Examine the tree associated with the model.
  hier_model.tree
  #Examine leaf cell types contained in the tree.
  hier_model.tree.cell_types(leaf_only = True)
  ```
  For details on what a tree is and how to construct, modify, and visualise it, see `4.`. For details on what a hierarchical model is and how to inspect, modify, and train it, see `5.`.
  </details>

+ <details>
  <summary><strong>1.4. Celltyping based on an LCPN hierarchical model</strong></summary>

  NN.
  </details>
</details>

<details>
<summary><strong>2. Use as the command line</strong></summary>
</details>

<details>
<summary><strong>3. Use as Docker/Singularity container</strong></summary>
</details>

<details>
<summary><strong>4. Construct, modify, and visualise a cell type hierarchy</strong></summary>

+ <details>
  <summary><strong>4.1. Structure of a cell type hierarchical tree</strong></summary>

  The basic unit of a cell type hierarchical tree is a node. In CellTypist, a node is represented by the [TreeNode](https://celltypist.readthedocs.io/en/latest/celltypist.tree.TreeNode.html) class, which contains information of a given cell type, including the mandatory attribute `original_name` (a unique cell type name), as well as optional attributes such as `cell_ontology_id`, `node_description`, `tissue_origin`, `markers`, `size`, `model`, and `children`. A JSON-like schematic looks like this:
  ```json
  {
    "original_name": "T cell",
    "cell_ontology_id": "CL:0000084",
    "node_description": "a type of lymphocyte responsible for cell-mediated immunity",
    "tissue_origin": ["blood", "lymphoid tissue"],
    "markers": ["CD3D"],
    "children": [
      {
        "original_name": "CD4+ T cell",
        "cell_ontology_id": "CL:0000624",
        "node_description": "helper T cell subtype",
        "markers": ["CD3D", "CD4"]
      },
      {
        "original_name": "CD8+ T cell",
        "cell_ontology_id": "CL:0000625",
        "node_description": "cytotoxic T cell subtype",
        "markers": ["CD3D", "CD8A"]
      }
    ]
  }
  ```
  Each node has a `children` attribute, which is either an empty list (`[]`) for leaf nodes or a list of [TreeNode](https://celltypist.readthedocs.io/en/latest/celltypist.tree.TreeNode.html) instances for internal nodes. In the above example, "CD4+ T cell" has no children (`children = []`), so the field is simply omitted in the JSON. Note that any child (e.g., "CD4+ T cell") may itself have children further, recursively extending the hierarchy.  

  Built upon [TreeNode](https://celltypist.readthedocs.io/en/latest/celltypist.tree.TreeNode.html), a cell type hierarchical tree in CellTypist is represented by the [Tree](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html) class: a `Tree` is a wrapper around a single `TreeNode` (the `root` node) together with a unique identifier (`handle`), plus optional fields for describing the tree.
  ```json
  {
    "handle": "T_Cell_Tree",
	"root": {
      "original_name": "T cell",
      "cell_ontology_id": "CL:0000084",
      "node_description": "a type of lymphocyte responsible for cell-mediated immunity",
      "tissue_origin": ["blood", "lymphoid tissue"],
      "markers": ["CD3D"],
      "children": [
        {
          "original_name": "CD4+ T cell",
          "cell_ontology_id": "CL:0000624",
          "node_description": "helper T cell subtype",
          "markers": ["CD3D", "CD4"]
        },
        {
          "original_name": "CD8+ T cell",
          "cell_ontology_id": "CL:0000625",
          "node_description": "cytotoxic T cell subtype",
          "markers": ["CD3D", "CD8A"]
        }
      ]
    }
  }
  ```
  Note that the value of `root` is exactly the "T cell" node shown earlier.
  </details>

+ <details>
  <summary><strong>4.2. Construct a cell type hierarchical tree</strong></summary>

  A simple way to construct a tree is to load it from an existing JSON file (see the second JSON schematic in `4.1.`) using [Tree.from_json](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.from_json).
  ```python
  from celltypist.tree import TreeNode, Tree
  tree = Tree.from_json('some_name.json')
  ```
  In principle, one could write the JSON file from scratch to capture the full hierarchy. However, for complex trees with deeply nested children, it is more practical to start from the root node (or ancestor clade) and gradually expand the structure by adding descendants. The example below shows the content of a starting JSON file (`root_only_tree.json`).
  ```json
  {
    "handle": "T_Cell_Tree",
	"root": {
      "original_name": "T cell",
      "node_description": "a type of lymphocyte responsible for cell-mediated immunity",
      "tissue_origin": ["blood", "lymphoid tissue"]
    }
  }
  ```
  In this file, you may omit any keys without values, except for the tree-level `handle` and `root`, and the node-level `original_name`. Next, load this file as a [Tree](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html) object.
  ```python
  #A tree with a single (root) node.
  tree = Tree.from_json('root_only_tree.json')
  ```
  An alternative way to initialise a `Tree` object, without relying on an external JSON file, is to use the [TreeNode](https://celltypist.readthedocs.io/en/latest/celltypist.tree.TreeNode.html) and [Tree](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html) constructors. This involves first creating a `TreeNode`, then initialising a `Tree`:
  ```python
  #Create a node.
  root_node = TreeNode(original_name = "T cell", node_description = "a type of lymphocyte responsible for cell-mediated immunity", tissue_origin = ["blood", "lymphoid tissue"])
  #Use this node as `root` to initialise a tree.
  tree = Tree(handle = "T_Cell_Tree", root = root_node)
  ```
  The resulting `tree` is equivalent to the one created with `Tree.from_json('root_only_tree.json')`. Note that in addition to `node_description` and `tissue_origin`, the `TreeNode` constructor also accepts other optional parameters such as `cell_ontology_id`, `markers`, `size`, and `model` (each explained in `4.3.`), as well as any custom fields. The `Tree` constructor likewise accepts custom fields besides the mandatory `handle` and `root`.  

  Next we create a "CD4+ T cell" node and a "CD8+ T cell" node, and add them as children of the root node "T cell".
  ```python
  #Create a "CD4+ T cell" node.
  CD4_node = TreeNode(original_name = "CD4+ T cell", cell_ontology_id = "CL:0000624", node_description = "helper T cell subtype", markers = ["CD3D", "CD4"])
  #Create a "CD8+ T cell" node.
  CD8_node = TreeNode(original_name = "CD8+ T cell", cell_ontology_id = "CL:0000625", node_description = "cytotoxic T cell subtype", markers = ["CD3D", "CD8A"])
  #Add both nodes as children of the root node "T cell".
  tree.add_children(CD4_node, CD8_node, parent = "T cell")
  ```
  Here the [add_children](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.add_children) method attaches the two new nodes as children of the designated parent node "T cell".  
  
  We can further create a new "γδ T cell" node and append it to the child list of "T cell" using the same method.
  ```python
  #Create a "γδ T cell" node.
  gamma_delta_node = TreeNode(original_name = "γδ T cell", node_description = "gamma-delta T cell subtype")
  #Append this node to the child list of "T cell".
  tree.add_children(gamma_delta_node, parent = "T cell")
  ```
  Now the root node "T cell" has three children. To preserve this structure, we can [write](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.write) the tree out as a JSON file so that our progress is not lost.
  ```python
  #Write out the tree locally.
  tree.write('T_cell_tree.json')
  ```
  This file can later be loaded back as the starting point for further tree extensions.
  ```python
  #Load the tree as a `Tree` instance.
  tree = Tree.from_json('T_cell_tree.json')
  ```
  Lastly, we create two nodes ("Naive CD4+ T cell" and "Memory CD4+ T cell") and add them as children of the "CD4+ T cell" node.
  ```python
  naive_cd4_node = TreeNode(original_name = "Naive CD4+ T cell", markers = ["CCR7", "SELL"])
  memory_cd4_node = TreeNode(original_name = "Memory CD4+ T cell")
  tree.add_children(naive_cd4_node, memory_cd4_node, parent = "CD4+ T cell")
  ```
  (Over)write the tree.
  ```python
  tree.write('T_cell_tree.json')
  ```
  </details>

+ <details>
  <summary><strong>4.3. Attributes and properties of a tree and its nodes</strong></summary>

  Load the tree built in `4.2.`.
  ```python
  tree = Tree.from_json('T_cell_tree.json')
  ```
  A schematic of its hierarchy looks like this:
  ```text
  T cell
  ├── CD4+ T cell
  │   ├── Naive CD4+ T cell
  │   └── Memory CD4+ T cell
  ├── CD8+ T cell
  └── γδ T cell
  ```
  This tree contains six nodes (four internal and two leaf nodes). You can retrieve any node, returned as a [TreeNode](https://celltypist.readthedocs.io/en/latest/celltypist.tree.TreeNode.html) instance, using the method [find_node](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.find_node).
  ```python
  #Retrieve the "CD4+ T cell" node.
  cd4_node = tree.find_node('CD4+ T cell')
  ```
  Each [TreeNode](https://celltypist.readthedocs.io/en/latest/celltypist.tree.TreeNode.html) has the following attributes, which can be accessed directly (e.g., `cd4_node.original_name`).
  <div align="center">

  |Attribute name  |Description                                                                                |Mandatory during [TreeNode](https://celltypist.readthedocs.io/en/latest/celltypist.tree.TreeNode.html) construction|Note     |
  |:---:           |:---:                                                                                      |:---:                                   |:---:                                                                               |
  |original_name   |The original and display name of the node (cell type)                                      |Yes                                     |Must be unique across all nodes within a tree                                       |
  |internal_name   |A programmatic version of the original name with special characters replaced by underscores|No                                      |Ignored for users                                                                   |
  |cell_ontology_id|Reference ID of the node (cell type) from the controlled Cell Ontology                     |No                                      |Empty string if not provided                                                        |
  |node_description|Description of the node (cell type)                                                        |No                                      |Empty string if not provided                                                        |
  |tissue_origin   |A list of tissue sources of the node (cell type)                                           |No                                      |Empty list if not provided                                                          |
  |markers         |A list of marker genes of the node (cell type)                                             |No                                      |Empty list if not provided                                                          |
  |size            |Number of cells contained in this node (cell type)                                         |No                                      |Typically populated during hierarchical model training. 0 if not provided           |
  |children        |A list of `TreeNode` instances representing *direct* child nodes (cell types) of the node  |No                                      |Typically added by [add_children](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.add_children) after initialisation. Empty list for leaf nodes|
  |model           |Path to a CellTypist model used for classifying child cell types of the given internal node|No                                      |Typically populated during hierarchical model training. Empty string if not provided|
  |any custom attr |N/A                                                                                        |No                                      |Do not overlap with above attributes                                                |
  </div>

  Note that the names of a parent node's immediate children can be accessed with `[child.original_name for child in cd4_node.children]`, or equivalently through the shortcut property `cd4_node.child_names`.  

  Built on [TreeNode](https://celltypist.readthedocs.io/en/latest/celltypist.tree.TreeNode.html), a [Tree](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html) defines three types of attributes:
  1. `handle` - a machine-friendly unique identifier for the tree
  2. `root` - a `TreeNode` serving as the root of the tree (N.B. a root is not necessarily a single terminal node, it can itself have recursive descendants)
  3. any custom attributes provided when using the [Tree](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html) constructor  

  Moreover, a `Tree` exposes several properties that can be accessed directly (e.g., `tree.depth`).
  <div align="center">

  |Property name   |Description                                              |Note                                 |
  |:---:           |:---:                                                    |:---:                                |
  |depth           |The depth of the tree                                    |The root has depth 1 (not 0)         |
  |n_leaves        |The number of leaf nodes contained in the tree           |Number of finest-grained cell types  |
  |n_nodes         |The number of total nodes contained in the tree          |Includes both leaf and internal nodes|
  |n_leaves_by_node|Dictionary mapping each node name to its total leaf count|Useful for reordering the tree       |
  </div>
  </details>

+ <details>
  <summary><strong>4.4. Query and manipulate a cell type hierarchical tree</strong></summary>

  Load the tree built in `4.2.`.
  ```python
  tree = Tree.from_json('T_cell_tree.json')
  ```
  A schematic of its hierarchy looks like this:
  ```text
  T cell
  ├── CD4+ T cell
  │   ├── Naive CD4+ T cell
  │   └── Memory CD4+ T cell
  ├── CD8+ T cell
  └── γδ T cell
  ```
  Although a `Tree` is just a wrapper around a `TreeNode`, it provides convenient methods for **querying** elements within the hierarchy:
  <div align="center">

  |Method name                                                                                                                                |Description                                       |Usage example                                              |Note                                                            |
  |:---:                                                                                                                                      |:---:                                             |:---:                                                      |:---:                                                           |
  |[cell_types](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.cell_types)                        |Return the cell type names contained in the tree  |`tree.cell_types(leaf_only = True)`                        |Setting `leaf_only = False` will return all cell type names     |
  |in                                                                                                                                         |Check whether a given cell type exists in the tree|`'CD4+ T cell' in tree`                                    |Check through all nodes in the tree, not just the leaf nodes    |
  |[iter_nodes](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.iter_nodes)                        |Iterate over nodes in the tree                    |`for node in tree.iter_nodes(leaf_only = False):`          |Setting `leaf_only = True` will iterate through only leaf nodes |
  |[find_node](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.find_node)                          |Find a node in the tree                           |`tree.find_node('CD4+ T cell')`                            |Return a `TreeNode` if found                                    |
  |[find_siblings](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.find_siblings)                  |Find a node's siblings which have the same parent |`tree.find_siblings('CD4+ T cell', return_names = True)`   |Setting `return_names = False` will return a list of `TreeNode`s|
  |[find_parent](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.find_parent)                      |Find the parent of a node in the tree             |`tree.find_parent('CD4+ T cell')`                          |Return `None` for the root node since it has no parent          |
  |[find_children](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.find_children)                  |Find the children of a node in the tree           |`tree.find_children('CD4+ T cell', return_names = True)`   |Setting `return_names = False` will return a list of `TreeNode`s|
  |[extract_subtree](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.extract_subtree)              |Extract a subtree rooted at a given node          |`tree.extract_subtree('CD4+ T cell')`                      |Return a new tree while leaving the original tree unchanged     |
  |[prune_by_depth](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.prune_by_depth)                |Prune the tree to a specified depth               |`tree.prune_by_depth(max_depth = 2)`                       |Return a new tree while leaving the original tree unchanged     |
  |[extract_path](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.extract_path)                    |Extract the path from the root to a given node    |`tree.extract_path('CD4+ T cell')`                         |Return a list of `TreeNode`s from the root to a given node      |
  |[lowest_common_ancestor](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.lowest_common_ancestor)|Find the lowest common ancestor of two nodes      |`tree.lowest_common_ancestor('CD4+ T cell', 'CD8+ T cell')`|Return the lowest common ancestor as a `TreeNode`               |
  </div>

  Note that the above methods do not alter the structure of the original tree.  

  The most important methods are arguably those for **manipulating** a tree, since a tree is continually subject to change and expansion. A full list of these methods is: 
  <div align="center">

  |Method name                                                                                                                    |Description                                                                    |Usage example                                                                                           |Note                                                                  |
  |:---:                                                                                                                          |:---:                                                                          |:---:                                                                                                   |:---:                                                                 |
  |[update](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.update)                    |Update attributes of a given node in the tree with provided keyword arguments  |`tree.update('CD4+ T cell', markers = ['CD3E', 'CD4'])`                                                 |Useful for modifying existing attributes or adding new attributes     |
  |[add_children](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.add_children)        |Add one or more nodes under a given parent node in the tree                    |`tree.add_children(new_node1, new_node2, parent = 'CD4+ T cell')`                                       |The tree is modified and the list of added nodes is returned          |
  |[add_node](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.add_node)                |Add a single node under a given parent node in the tree                        |`tree.add_node(new_node, parent = 'CD4+ T cell')`                                                       |Similar to `tree.add_children(new_node, parent = 'CD4+ T cell')`      |
  |[remove_node](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.remove_node)          |Remove a node from the tree                                                    |`tree.remove_node('CD4+ T cell')`                                                                       |The tree is modified and the removed node is returned                 |
  |[move_node](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.move_node)              |Move a node (and its descendants) from its current parent to a new parent      |`tree.move_node('Naive CD4+ T cell', to = 'CD8+ T cell')`                                               |The tree is modified and the moved node is returned                   |
  |[replace_node](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.replace_node)        |Replace a node in the tree (and its descendants) with a new node               |`tree.replace_node('CD4+ T cell', by = new_node)`                                                       |The tree is modified and the replaced (old) node is returned          |
  |[reorder_children](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.reorder_children)|Reorder the direct children of a parent node to match a new order              |`tree.reorder_children(parent = 'CD4+ T cell', new_order = ['Memory CD4+ T cell', 'Naive CD4+ T cell'])`|The tree is modified and the reordered list of nodes is returned      |
  |[sort_children](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.sort_children)      |Sort the direct children of a parent node by the number of leaves they contain |`tree.sort_children(parent = 'CD4+ T cell', recursive = True, descending = True)`                       |The tree is modified and the sorted list of nodes is returned         |
  |[sort_tree](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.sort_tree)              |Sort the children of each node in the tree by the number of leaves they contain|`tree.sort_tree(recursive = True, descending = True)`                                                   |Setting `recursive = False` will sort by the number of direct children|
  |[assign_size](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.assign_size)          |Assign cell counts to all nodes based on the provided leaf-level annotations   |`tree.assign_size(leaf_anno = a_leaf_label_vector)`                                                     |Equivalent to counting cells contained in the subtree of each node    |
  |[copy](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.copy)                        |Create a deep copy of the tree                                                 |`copied_tree = tree.copy()`                                                                             |Changes made to `copied_tree` will not impact the original `tree`     |
  </div>

  Note that except for `copy`, the above methods modify the structure of the original tree in place. Given this, it is good practice to run the [validate](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.validate) method *regularly* to sanity-check the entire tree (e.g., ensuring updated attributes have the correct type and names of newly added nodes do not duplicate existing ones in the tree). This is particularly useful as a health check before saving the tree as a JSON file.
  ```python
  #Sanity-check the entire tree. Raises a corresponding error if found.
  tree.validate()
  #Write out the tree locally.
  tree.write('T_cell_tree.json')
  ```
  </details>

+ <details>
  <summary><strong>4.5. Generate multi-level cell type annotations from leaf labels</strong></summary>

  In single-cell datasets, the cell metadata table often contains a column of fine-grained annotations (i.e., leaf-level labels). By coupling these labels with a predefined cell type hierarchy, they can be mapped to broader categories at different levels. For example, consider the tree below:
  ```text
			   A
			   │
		┌──────┼──────┐
		│      │      │
	  A1     A2     A3
             		  │
             	   ┌──┴──┐
             	   │     │
             	  A3a   A3b
  ```
  This tree has a total depth of three, corresponding to three levels of cell type annotations. At level 3 (the leaf level), if five cells are annotated as `A1 A3a A3a A2 A3b`, their corresponding level 2 and level 1 annotations are as follows:
  ```text
  Level 3: A1  A3a  A3a  A2  A3b
  Level 2: A1  A3   A3   A2  A3
  Level 1: A   A    A    A   A
  ```
  At level 2, all level 3 annotations are converted into their level 2 ancestors. This process continues until all three levels of annotations are obtained.
  ### Why do this?
  > This multi-level mapping is a prerequisite for hierarchical training (detailed in `5.`). As a preview, we can either 1) train a global level-2 classifier based on the level-2 annotation vector, or 2) train a local classifier to distinguish A3a vs. A3b cells within A3 (using the level-2 vector to locate A3 cells and the level-3 vector to locate A3a/A3b cells).  

  CellTypist provides this functionality through the [get_multilevel_anno](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.get_multilevel_anno) method, which takes a vector of leaf-level annotations (e.g., a `list`, `tuple`, or `Series`) as input. Using the example `tree` above:
  ```python
  leaf_anno = ["A1", "A3a", "A3a", "A2", "A3b"]
  multi_anno = tree.get_multilevel_anno(leaf_anno)
  ```
  The output (`multi_anno`) is a `DataFrame` with one column per annotation level.
  <div align="center">

  | level_1_anno | level_2_anno | level_3_anno |
  |:------------:|:------------:|:------------:|
  | A            | A1           | A1           |
  | A            | A3           | A3a          |
  | A            | A3           | A3a          |
  | A            | A2           | A2           |
  | A            | A3           | A3b          |
  </div>

  Note that if the input is a `Series`, this output `DataFrame` preserves its index (i.e., cell names). Also, there is a `prefix` parameter for [get_multilevel_anno](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.get_multilevel_anno) by which you can add a custom string to the beginning of each column name.
  </details>

+ <details>
  <summary><strong>4.6. Visualise a cell type hierarchical tree</strong></summary>

  A [Tree](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html) can be visualised using the [celltypist.treeviz](https://celltypist.readthedocs.io/en/latest/celltypist.treeviz.html) function:
  ```python
  #Visualise a tree directly.
  celltypist.treeviz(tree)
  #Or alternatively, save the tree plot.
  celltypist.treeviz(tree, show = False, save = 'tree.pdf')
  ```
  By default, this plot uses a diagonal layout (`layout = "diagonal"`). Setting `layout = "rectangular"` will instead draw the tree with a rectangular layout.
  ```text
  layout = "diagonal"      layout = "rectangular"
           ●                        ───●
          ╱                        |
         ●───●                     ●───●
          ╲                        |
           ●                        ───●
  ```
  With the default `direction = "right"`, the tree grows towards the right. Setting `direction = "down"` will make the tree grow downwards.
  ```text
  direction = "right"      direction = "down"
           ●       
          ╱                     ───●───
         ●───●                 |   |   |
          ╲                    ●   ●   ●
           ●                      
  ```
  Since the order of children for a given parent node can be rearranged without affecting the tree, you can freely adjust the order of children (for example, using the method [reorder_children](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html#celltypist.tree.Tree.reorder_children) as described in `4.4.`). To quickly organise the entire tree, one option is to order the children by their complexity. Specifically, when `sort = True` (the default is `sort = False` as we expect an already-ordered tree in most cases), each internal node’s children are arranged from most to least complex.
  ```python
  #Reorder the children of every internal node and visualise the tree.
  celltypist.treeviz(tree, sort = True)
  #Node complexity = total descendant leaves (`recursive = True`, the default) or number of direct children (`recursive = False`).
  celltypist.treeviz(tree, sort = True, recursive = False)
  #Reorder children of each internal node from least to most complex with `descending = False` (default to `descending = True`).
  celltypist.treeviz(tree, sort = True, descending = False)
  ```
  You can also choose whether to display cell type labels for internal nodes (default: `show_node_label = False`) and leaf nodes (default: `show_leaf_label = True`). Other parameters for controlling the shapes, colors, sizes, and alignments of branches, nodes, and labels can be found in [celltypist.treeviz](https://celltypist.readthedocs.io/en/latest/celltypist.treeviz.html).
  </details>
</details>

<details>
<summary><strong>5. Inspect, modify, and train a hierarchical model</strong></summary>

+ <details>
  <summary><strong>5.1. Structure of a hierarchical model</strong></summary>

  The hierarchical model in CellTypist is implemented through the [HierModel](https://celltypist.readthedocs.io/en/latest/celltypist.models.HierModel.html) class.  

  Load the default hierarchical model (`Human_Tissue_Immune_LCPN.pkl`), which includes various immune and hematopoietic cell types across human tissues.
  ```python
  #If the `model` argument is not provided, it will default to `Human_Tissue_Immune_LCPN.pkl`.
  hier_model = HierModel.load()
  #Show summary information of the model.
  hier_model
  ```
  Each [HierModel](https://celltypist.readthedocs.io/en/latest/celltypist.models.HierModel.html) has an associated [Tree](https://celltypist.readthedocs.io/en/latest/celltypist.tree.Tree.html), in which all cell types are represented as either internal or leaf nodes.
  ```python
  #Access the tree of the model.
  hier_model.tree
  #You can visualise the tree using `celltypist.treeviz`, as detailed in `4.6.`.
  celltypist.treeviz(tree)
  ```
  A hierarchical model is an ensemble model composed of multiple flat models, operating in one of two modes: Local Classifier per Parent Node (`LCPN`), which trains a local classifier for each parent node (coarse cell type), and Local Classifier per Level (`LCL`), which trains a local classifier at each level (depth) of the tree.
  ```python
  #Examine the mode in which the model is trained.
  hier_model.mode
  ```
  In the LCPN mode, each parent node (a coarse-level cell type) often contains a local classifier trained to distinguish among its child cell types (subtypes). These classifiers are stored as flat CellTypist models, organised in the attribute `hier_model.model_mapping` as a dictionary.
  ```python
  #Access all local classifiers of an LCPN hierarchical model.
  hier_model.model_mapping
  ```
  Specifically, each parent node (which is a [TreeNode](https://celltypist.readthedocs.io/en/latest/celltypist.tree.TreeNode.html) instance) points to its classifier through the attribute `model` (`node.model`), and `hier_model.model_mapping` acts as the lookup table that stores all classifier objects. For example, to retrieve the local classifier linked to `T cell`, we first locate the `T cell` node in the tree.
  ```python
  #Fetch the 'T cell' tree node.
  T_node = hier_model.tree.find_node('T cell')
  ```
  Examine the filename/key of `T cell`'s associated classifier.
  ```python
  model_key = T_node.model
  print(model_key)
  # -> Output: 'T_cell.pkl'
  ```
  Retrieve the CellTypist flat model (i.e., local classifier) linked to this node.
  ```python
  T_local_model = hier_model.model_mapping[model_key]
  ```
  This local classifier is intended for classifying the *direct* child cell types of `T cell`. To check which cell types these are, access `T_local_model.cell_types` or `T_node.child_names`.
  ```python
  #Shows the subtypes of `T cell` in the tree.
  T_local_model.cell_types
  #Alternatively
  T_node.child_names
  ```
  In short, the node’s `model` attribute stores the name of the classifier, while `hier_model.model_mapping` keeps the actual [Model](https://celltypist.readthedocs.io/en/latest/celltypist.models.Model.html) object by mapping the classifier name.  
  </details>
</details>

<details>
<summary><strong>Supplemental guidance</strong></summary>
</details>

# Citation
Dominguez Conde et al., Cross-tissue immune cell analysis reveals tissue-specific features in humans. Science 376, eabl5197 (2022). [Link](https://doi.org/10.1126/science.abl5197)
