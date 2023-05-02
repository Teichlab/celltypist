<p align="left"><img src="https://github.com/Teichlab/celltypist/blob/main/docs/source/_static/img/logo_celltypist.png" width="250" height="85"></p>

[![Python Versions](https://img.shields.io/badge/python-3.6+-brightgreen.svg)](https://pypi.org/project/celltypist) [![Documentation Status](https://readthedocs.org/projects/celltypist/badge/?version=latest)](https://celltypist.readthedocs.io/en/latest/?badge=latest)

CellTypist is an automated tool for cell type classification, harmonisation, and integration.
- _classification_: transfer cell type labels from the reference to query dataset
- _harmonisation_: match and harmonise cell types defined by independent datasets
- _integration_: integrate cell and cell types with supervision from harmonisation

# CellTypist website
Information of CellTypist can be also found in our CellTypist portal. [![Website www.celltypist.org](https://img.shields.io/website-up-down-brightgreen-red/http/shields.io.svg)](https://www.celltypist.org/)

# Interactive tutorials
### Classification
[Using CellTypist for cell type classification ![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Teichlab/celltypist/blob/main/docs/notebook/celltypist_tutorial.ipynb)  
[Using CellTypist for multi-label classification ![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Teichlab/celltypist/blob/main/docs/notebook/celltypist_tutorial_ml.ipynb)  
[Best practice in large-scale cross-dataset label transfer using CellTypist ![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Teichlab/celltypist/blob/main/docs/notebook/celltypist_tutorial_cv.ipynb)
### Harmonisation
[Using CellTypist for cell type harmonisation ![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Teichlab/celltypist/blob/main/docs/notebook/celltypist_tutorial_harmonisation.ipynb)
### Integration
[Using CellTypist for annotation-aware data integration ![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Teichlab/celltypist/blob/main/docs/notebook/celltypist_tutorial_integration.ipynb)

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
</details>

# Usage (harmonisation)

<details>
<summary><strong>1. Cross-dataset cell type harmonisation</strong></summary>

+ <details>
  <summary><strong>1.1. Cell type harmonisation</strong></summary>

  The input [AnnData](https://anndata.readthedocs.io/en/latest/) needs two columns in `.obs` representing dataset origin and cell original annotation respectively. The aim is to harmonise cell types across datasets using [celltypist.harmonize](https://celltypist.readthedocs.io/en/latest/celltypist.harmonize.html).  
    
  Internally, transcriptional distances between cells and cell types (denoted here as the cell centroid) will first be calculated. Since cell type is usually defined at the cluster level and no cluster is 100% pure, you can set `filter_cells = True` (default to `False`) to filter out cells whose gene expression profiles do not correlate most with the cell type they belong to. This will speed up the run as only a subset of cells are used, but will render these filtered cells unannotated (see `2.2.`). Distances are calculated at either gene or low-dimensional space. The latter is preferred to denoise the data by providing a latent representation via the argument `use_rep` (default to PCA coordinates).
  ```python
  #`use_rep` can be omitted here as it defaults to 'X_pca'.
  alignment = celltypist.harmonize(adata, dataset = 'dataset_column', cell_type = 'celltype_column', use_rep = 'X_pca')
  ```
  If `X_pca` is not detected in `.obsm` and no other latent representations are provided via `use_rep`, gene expression matrix in `.X` will be used to calculate the distances. In such case, subsetting the AnnData to informative genes (e.g. highly variable genes) is suggested and `.X` should be log-normalised (to a constant total count per cell).  
    
  The resulting `alignment` is an instance of the class [DistanceAlignment](https://celltypist.readthedocs.io/en/latest/celltypist.contro.align.DistanceAlignment.html) as defined by CellTypist, and can be written out as follows.
  ```python
  #Save the harmonisation output.
  alignment.write('/path/to/local/folder/some_name.pkl')
  ```
  </details>

+ <details>
  <summary><strong>1.2. Cell type harmonisation with PCT</strong></summary>

  Inferring cell type relationships based on directly calculated distances will suffice in most cases due to a normalisation procedure applied to the derived distances. If a very strong batch effect exists across datasets, you can turn on `use_pct = True` (default to `False`) to predict instead of calculate these distances. Through this parameter, a predictive clustering tree (PCT) is built for each dataset, and distances between cells in query datasets and cell types in the reference dataset are predicted, often resulting in unbiased distance measures.
  ```python
  #Use PCT to predict transcriptional cell-cell distances across datasets.
  alignment = celltypist.harmonize(adata, dataset = 'dataset_column', cell_type = 'celltype_column', use_rep = 'X_pca', use_pct = True)
  ```
  Due to the nonparametric nature of PCT, the format of the expression `.X` in the AnnData is flexible (normalised, log-normalised, z-scaled, etc.), but subsetting the AnnData to highly variable genes is always suggested. To avoid overfitting, each PCT is pruned at nodes where no further splits are needed based on F-test, which is turned on by default (`F_test_prune = True`). You can increase the p-value cutoff (default to 0.05, `p_thres = 0.05`) to prune fewer nodes for improved accuracy at the cost of reduced generalisability.
  </details>

+ <details>
  <summary><strong>1.3. Specify the dataset order</strong></summary>

  In CellTypist, datasets are iteratively incorporated and harmonised. The order of datasets can be specified by providing a list of dataset names to the argument `dataset_order`. Otherwise, the order will be determined by CellTypist through iteratively adding a dataset that is most similar (i.e., more shared cell types) to the datasets already incorporated. This behaviour can be disabled by setting `reorder_dataset = False` (default to `True`) and an alphabetical order of datasets will be used.
  ```python
  #Specify the order of datasets to be harmonised.
  alignment = celltypist.harmonize(adata, dataset = 'dataset_column', cell_type = 'celltype_column', use_rep = 'X_pca', dataset_order = a_list_of_datasets)
  ```
  </details>

+ <details>
  <summary><strong>1.4. Categories of harmonised cell types</strong></summary>

  Four kinds of harmonisations are anchored with [celltypist.harmonize](https://celltypist.readthedocs.io/en/latest/celltypist.harmonize.html):
     1) Novel cell types as determined by `maximum_novel_percent` (default to `0.05`). In each harmonisation iteration, a cell type (or meta-cell-type) whose maximal alignment fraction is < `maximum_novel_percent` with any cell types in any other datasets is designated as a novel cell type (`NONE`).
     2) One-to-one aligned cell types as determined by `minimum_unique_percents` and `minimum_divide_percents`. If the alignments (in both directions) between two cell types from two respective datasets are greater than `minimum_unique_percents`, plus that these alignments are not one-to-many (see the third point below), this will be an 1:1 (`=`) match. Dynamic thresholds of `minimum_unique_percents` (default to 0.4, 0.5, 0.6, 0.7, 0.8) and `minimum_divide_percents` (default to 0.1, 0.15, 0.2) are exhaustively tested until the least number of alignments is found between datasets.
     3) One-to-many (or many-to-one) aligned cell types as determined by `minimum_unique_percents` and `minimum_divide_percents`. If one cell type has more than two cell types aligned in the other dataset with a match proportion greater than `minimum_divide_percents`, and these matched cell types have a back-match proportion greater than `minimum_unique_percents`, this will be an 1:N (`∋`) or N:1 (`∈`) match. Dynamic thresholds of `minimum_unique_percents` (default to 0.4, 0.5, 0.6, 0.7, 0.8) and `minimum_divide_percents` (default to 0.1, 0.15, 0.2) are exhaustively tested until the least number of alignments is found between datasets.
     4) Unharmonised cell types. If after the above categorisation, a cell type remains unharmonised, then this cell type will be an unharmonised cell type (`UNRESOLVED`).  
    
  |If there are many datasets to harmonise and each dataset has many cell types, harmonisation may take longer time. You can restrict the test scope of `minimum_unique_percents` and `minimum_divide_percents` to reduce runtime. The default is a 15 (5X3) combo test; setting the two parameters to, for example a 3X2 combo, can decrease 60% of the runtime.|
  |:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
  ```python
  #`minimum_unique_percents` is set to three values (default is 0.4, 0.5, 0.6, 0.7, 0.8).
  #`minimum_divide_percents` is set to two values (default is 0.1, 0.15, 0.2).
  alignment = celltypist.harmonize(adata, dataset = 'dataset_column', cell_type = 'celltype_column', use_rep = 'X_pca', minimum_unique_percents = [0.5, 0.6, 0.7], minimum_divide_percents = [0.1, 0.15])
  ```
  </details>
</details>

<details>
<summary><strong>2. Inspection of the harmonisation result</strong></summary>

+ <details>
  <summary><strong>2.1. Harmonisation table</strong></summary>

  The previously saved harmonisation object can be loaded using `celltypist.DistanceAlignment.load`.
  ```python
  alignment = celltypist.DistanceAlignment.load('/path/to/local/folder/some_name.pkl')
  ```
  In `alignment`, the harmonisation table, which summarises cell types across datasets into semantically connected ones, is stored as the attribute `.relation` (`alignment.relation`). One illustrative example is:
  <div align="center">

  |D1   |relation|D2   |relation|D3        |
  |:---:|:---:   |:---:|:---:   |:---:     |
  |A    |=       |B    |=       |C         |
  |D    |=       |NONE |=       |UNRESOLVED|
  |E    |∈       |G    |=       |H         |
  |F    |∈       |G    |=       |I         |
  |J    |=       |K    |∋       |L         |
  |J    |=       |K    |∋       |M         |
  </div>

  The table columns are the dataset1 name, relation, dataset2 name, ..., all the way to the name of the last dataset. Accordingly, each row of the table is a list of cell types connected by predefined symbols of `=`, `∈`, and `∋`. In addition to cell type names, there are two extra definitions of `NONE` and `UNRESOLVED` in the table, representing two levels of novelties (see `1.4.`).  
    
  The table should be interpreted from left to right. For example, for the first row `A = B = C`, although it may look like an 1:1 match between A and B plus an 1:1 match between B and C, a correct interpretation should be an 1:1 match between A and B, resulting in a meta cell type of `A = B`. This meta cell type, as a whole, has an 1:1 match with C, further leading to `A = B = C`. Similarly, for the second row `D = NONE = UNRESOLVED`, instead of a novel cell type D in dataset1, this cell type should be read as a dataset1-specific cell type not existing in dataset2 (`D = NONE`), which as a whole is unharmonised when aligning with dataset3 (`D = NONE = UNRESOLVED`).  
    
  Extending this interpretation to the third and fourth rows, they denote two cell types (E and F) in dataset1 collectively constituting the cell type G in dataset2. The resulting subtypes (`E ∈ G` and `F ∈ G`) are 1:1 matched with H and I in dataset3, respectively. For the last two rows, they describe the subdivision of a meta cell type (`J = K`) into L and M in dataset3, being more than a subdivision of K.  
    
  In the table, each row corresponds to a harmonised low-hierarchy cell type, in other words, the most fine-grained level of annotation that can be achieved by automatic alignment. At a high hierarchy, some cell types such as `E ∈ G = H` and `F ∈ G = I` belong to the same group. CellTypist defines a high-hierarchy cell type as fully connected rows in the harmonisation table. As a result, each high-hierarchy cell type is a cell type group independent of each other. This information can be accessed in the attribute `.groups` which is an array/vector with an length of the number of rows in the harmonisation table.
  ```python
  #Access the high-hierarchy cell types (cell type groups).
  alignment.groups
  ```
  </details>

+ <details>
  <summary><strong>2.2. Cell reannotation</strong></summary>

  After cell type harmonisation, each cell can be assigned a cell type label corresponding to a given row of the harmonisation table, denoted as the process of cell reannotation. By default, reannotation is enabled (`reannotate = True`) when using [celltypist.harmonize](https://celltypist.readthedocs.io/en/latest/celltypist.harmonize.html) and information of reannotated cell types is already in place as the attribute `.reannotation`.
  ```python
  #Access the cell reannotation information.
  alignment.reannotation
  ```
  This is a data frame with an example shown below. Unless `filter_cells = True` is set (see `1.1.`), all cells in the AnnData will be present in this data frame.
  <div align="center">

  |     |dataset|cell_type|reannotation         |group |
  |:---:|:---:  |:---:    |:---:                |:---: |
  |cell1|D1     |A        |A = B = C            |Group1|
  |cell2|D1     |D        |D = NONE = UNRESOLVED|Group2|
  |cell3|D2     |G        |E ∈ G = H            |Group3|
  |cell4|D2     |G        |F ∈ G = I            |Group3|
  |cell5|D3     |L        |J = K ∋ L            |Group4|
  |cell6|D3     |M        |J = K ∋ M            |Group4|
  </div>

  The four columns represent information of dataset origin, original author annotation, reannotated low- and high-hierarchy annotation, respectively. For the last column, it contains grouping (high-hierarchy) information, and each group corresponds to a subset of the harmonisation table. You can check this correspondence by coupling the table (`alignment.relation`) with the grouping (`alignment.groups`) (see `2.1.`).  
    
  Of note, due to several reasons including the clustering impurity of homogeneous cell populations, not all cells can be reannotated into low-hierarchy cell types, leading to `UNASSIGNED` cells within the `reannotation` column. The `group` column, however, is always populated with meaningful groups as this column is directly based on the harmonisation table.
  </details>

+ <details>
  <summary><strong>2.3. Meta-analysis</strong></summary>

  A distance matrix-like instance, which is from the class [Distance](https://celltypist.readthedocs.io/en/latest/celltypist.contro.distance.Distance.html) as defined by CellTypist, is also stashed in `alignment` as the attribute `.base_distance`.
  ```python
  #Access the distance object.
  alignment.base_distance
  ```
  The main content of this object is the distance matrix (`alignment.base_distance.dist_mat`) between all cells (rows) and all cell types (columns). Values in this matrix are either calculated (the default) or inferred (if `use_pct` is `True`) by `celltypist.harmonize`, and after a normalisation procedure, lie between 0 and 1. If there are strong cross-dataset batches, an inferred distance matrix obtained from the PCT algorithm is usually more accurate. Metadata of cells and cell types for this matrix can be found in `alignment.base_distance.cell` and `alignment.base_distance.cell_type`, which record raw information such as the dataset origin and original author annotation.  
    
  During the internal harmonisation process, each cell is assigned the most similar cell type from each dataset. This result is stored in the assignment matrix (`alignment.base_distance.assignment`), with rows being cells (cell metadata can be found in `alignment.base_distance.cell` as mentioned above), columns being datasets, and elements being the assigned cell types in different datasets. This matrix can be interpreted as a summary of multi-data label transfers.
  ```python
  #Access the cell type assignment result.
  alignment.base_distance.assignment
  ```
  Each column (corresponding to one dataset) of the assignment matrix can be thought as a unified naming schema when all cells are named by this given dataset.  
    
  CellTypist provides a quick way to summarise the above information including cells' distances and assignments into meta-analysis at the cell type level. Specifically, a distance matrix among all cell types can be obtained by:
  ```python
  #Get the cell-type-to-cell-type distance matrix.
  alignment.base_distance.to_meta()
  ```
  An optional `turn_binary = True` (default to `False`) can be added to turn the distance matrix into a cell membership matrix before meta-analysis, showing how cell types are assigned across datasets.
  ```python
  #Get the cell-type-to-cell-type membership matrix.
  alignment.base_distance.to_meta(turn_binary = True)
  ```
  </details>
</details>

<details>
<summary><strong>3. Reharmonisation</strong></summary>

+ <details>
  <summary><strong>3.1. Change the dataset order</strong></summary>

  The order of datasets used by `celltypist.harmonize` can be found in the attribute `.dataset_order` (`alignment.dataset_order`), which is either auto-determined by CellTypist or specified by the user (via the `dataset_order` parameter in `celltypist.harmonize`). This order is also reflected by the column order of the harmonisation table.  
    
  Along the order of datasets, optimal choices of `minimum_unique_percents` and `minimum_divide_percents` (see `1.4.`) in each iteration can be found in `alignment.minimum_unique_percents` and `alignment.minimum_divide_percents`. For instance, harmonising five datasets requires four iterations, and thus both `.minimum_unique_percents` and `.minimum_divide_percents` have a length of four.  
    
  CellTypist provides a method [best_align](https://celltypist.readthedocs.io/en/latest/celltypist.contro.align.DistanceAlignment.html#celltypist.contro.align.DistanceAlignment.best_align) to change the order of datasets post-harmonisation. Through this, datasets will be reharmonised in a different order (this post-harmonisation adjustment is more efficient than re-running `celltypist.harmonize` with a new order).
  ```python
  #Reharmonise cell types across datasets with a different dataset order.
  alignment.best_align(dataset_order = a_list_of_new_dataset_order)
  ```
  As in `celltypist.harmonize`, the combos of `minimum_unique_percents` and `minimum_divide_percents` will be tested to find the best alignment in each iteration. Importantly, as well as a full dataset list, you can provide a subset of datasets for reharmonisation. This is useful in terms of focusing on part of the data for inspection or visualisation (see `4.`).
  ```python
  #Reharmonise cell types across datasets with part of datasets.
  alignment.best_align(dataset_order = a_subset_of_dataset_names)
  ```
  A new harmonisation table will be generated in `alignment.relation`, which only includes datasets specified in `.best_align`. `.minimum_unique_percents` and `.minimum_divide_percents` are also overridden by new values used during reharmonisation.
  </details>

+ <details>
  <summary><strong>3.2. Reannotation</strong></summary>

  After changing the dataset order and reharmonising cell types, cells need to be reannotated based on the newly generated harmonisation table using the method [reannotate](https://celltypist.readthedocs.io/en/latest/celltypist.contro.align.DistanceAlignment.html#celltypist.contro.align.DistanceAlignment.reannotate).
  ```python
  #Reannotate cells based on the new harmonisation table.
  alignment.reannotate()
  ```
  Similarly, information of reannotated cells is stored in `alignment.reannotation`.
  </details>
</details>

<details>
<summary><strong>4. Visualisation</strong></summary>

+ <details>
  <summary><strong>4.1. Tree plot</strong></summary>

  The most intuitive way to visualise the harmonised cell types is the tree plot using the function [celltypist.treeplot](https://celltypist.readthedocs.io/en/latest/celltypist.treeplot.html).
  ```python
  #Visualise the harmonisation result with a tree plot.
  celltypist.treeplot(alignment)
  ```
  Alternatively, since only the harmonisation table (`alignment.relation`) is used when plotting this tree, `celltypist.treeplot` also accepts the input directly from the table. This is more convenient as a table is easier to manipulate, such as writing it out as a csv file and loading it later for tree plot.
  ```python
  #Write out the harmonisation table as a csv file.
  #Note - if cell type names contain commas, set a different `sep` here.
  alignment.relation.to_csv('/path/to/local/folder/HT.csv', sep = ',', index = False)
  ```
  ```python
  #Read the harmonisation table.
  HT = pd.read_csv('/path/to/local/folder/HT.csv', sep = ',')
  #Visualise the harmonisation result with a tree plot.
  celltypist.treeplot(HT)
  #Visualise the harmonisation result only for cell types (rows) of interest.
  celltypist.treeplot(HT[row_flag])
  ```
  In a tree plot, each column is a dataset and cell types are connected across datasets. By default, cell types belonging to one low hierarchy (one row in the harmonisation table) are in the same color. You can change the color scheme by providing a data frame to the `node_color` parameter, with three consecutive columns representing dataset, cell type, and color (in hex code), respectively. `node_color` can also be a data frame with columns of dataset, cell type, and numeric value (for mapping color gradient in combination with `cmap`). Other parameters controlling the appearance of the tree plot (node shape, line width, label size, figure size, etc.) are detailed in [celltypist.treeplot](https://celltypist.readthedocs.io/en/latest/celltypist.treeplot.html).  
  |The tree plot considers all pairs of reference-to-query assignments. Therefore, a restricted representation in two dimensionalities may overlay some cell types when they have complex 1:1 and 1:N intersections. These cross-connections are usually not solvable at 2D space; you may need to revisit the harmonisation table in some cases.|
  |:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
    
  By changing the dataset (column) order in each high-hierarchy cell type, broader (more divisible) cell types can be positioned to the left, followed by fine-grained cell types to the right. The resulting plot shows how different authors group these cell types, thereby being more characteristic of the potential underlying biological hierarchy. This hierarchy can be generated and visualised by adding `order_dataset = True`.
  ```python
  #Visualise the cell type hierarchy.
  #Again, the input can also be a harmonisation table.
  celltypist.treeplot(alignment, order_dataset = True)
  ```
  Because each high-hierarchy cell type is independent of each other, the new orders of datasets will be different across groups. To recognise the dataset origin of each cell type within the hierarchy, you can assign the same color or shape to cell types from the same dataset using the parameter `node_color` or `node_shape`. An example is:
  ```python
  #Cell types from the same dataset are in the same shape.
  #`node_shape` should be the same length as no. datasets in the harmonisation table.
  celltypist.treeplot(alignment, order_dataset = True, node_shape = list_of_shapes)
  ```
  Export the plot if needed.
  ```python
  celltypist.treeplot(alignment, show = False, save = '/path/to/local/folder/some_name.pdf')
  ```
  </details>

+ <details>
  <summary><strong>4.2. Sankey plot</strong></summary>

  The other way to visualise harmonised cell types is the Sankey plot by [celltypist.sankeyplot](https://celltypist.readthedocs.io/en/latest/celltypist.sankeyplot.html). CellTypist builds this plot on the [plotly](https://pypi.org/project/plotly) package. `plotly` is not mandatory when installing CellTypist, so you need to install it first if you want a visualisation form of Sankey diagram (and engines for exporting images such as [kaleido](https://pypi.org/project/kaleido)).
  ```python
  #Visualise the harmonisation result with a Sankey plot.
  #As with the tree plot, the input can also be a harmonisation table.
  celltypist.sankeyplot(alignment)
  ```
  Similar to the tree plot, this diagram shows how cell types are connected across datasets. Parameters controlling the appearance of the Sankey plot (node color, link color, figure size, etc.) are detailed in [celltypist.sankeyplot](https://celltypist.readthedocs.io/en/latest/celltypist.sankeyplot.html).  
    
  Different from the tree plot where novel (`NONE`) and unharmonised (`UNRESOLVED`) cell types are blank, in the Sankey plot they are colored in white and light grey, respectively. You can adjust these by changing the values of `novel_node_color` and `remain_node_color`.  
    
  Export the plot if needed.
  ```python
  #Export the image into html.
  celltypist.sankeyplot(alignment, show = False, save = '/path/to/local/folder/some_name.html')
  #Export the image into pdf.
  celltypist.sankeyplot(alignment, show = False, save = '/path/to/local/folder/some_name.pdf')
  ```
  </details>
</details>

# Usage (integration)

<details>
<summary><strong>1. Supervised data integration</strong></summary>

+ <details>
  <summary><strong>1.1. Specify batch and biological covariates</strong></summary>

  The input [AnnData](https://anndata.readthedocs.io/en/latest/) needs two columns in `.obs` representing the batch confounder and unified cell annotation respectively. The aim is to integrate cells by correcting batches and preserving biology (cell annotation) using [celltypist.integrate](https://celltypist.readthedocs.io/en/latest/celltypist.integrate.html).
  ```python
  #Integrate cells with `celltypist.integrate`.
  celltypist.integrate(adata, batch = 'a_batch_key', cell_type = 'a_celltype_key')
  ```
  With this function, CellTypist will build the neighborhood graph by searching neighbors across matched cell groups in different batches, on the basis of a low-dimensional representation provided via the argument `use_rep` (default to PCA coordinates).
  ```python
  #`use_rep` can be omitted here as it defaults to 'X_pca'.
  celltypist.integrate(adata, batch = 'a_batch_key', cell_type = 'a_celltype_key', use_rep = 'X_pca')
  ```
  The batch confounder can be the dataset origin, donor ID, or any relevant covariate. For the biological factor, it is the consistent annotation across cells, such as manual annotations of all cells, transferred cell type labels from a single reference model, and as an example here, the harmonised cell types from the CellTypist harmonisation pipeline (see the harmonisation section). Specifically, you can add two extra columns in the `.obs` of the input AnnData using the reannotation information from `alignment.reannotation`.
  ```python
  #Insert low- and high-hierarchy annotations into the AnnData.
  adata.obs[['harmonized_low', 'harmonized_high']] = alignment.reannotation.loc[adata.obs_names, ['reannotation', 'group']]
  ```
  Perform data integration using either of the two annotation columns.
  ```python
  #Integrate cells using the reannotated high-hierarchy cell annotation.
  celltypist.integrate(adata, batch = 'a_batch_key', cell_type = 'harmonized_high')
  #Not run; integrate cells using the reannotated low-hierarchy cell annotation.
  #celltypist.integrate(adata, batch = 'a_batch_key', cell_type = 'harmonized_low')
  ```
  Finally, generate a UMAP based on the reconstructed neighborhood graph.
  ```python
  sc.tl.umap(adata)
  ```
  </details>

+ <details>
  <summary><strong>1.2. Adjust the influence of annotation on integration</strong></summary>

  Influence of cell annotation on the data structure can range from forcibly merging the same cell types to a more lenient cell grouping. This is achieved by adjusting the parameter `n_meta_neighbors`.
  ```python
  #Actually the default value of `n_meta_neighbors` is 3.
  celltypist.integrate(adata, batch = 'a_batch_key', cell_type = 'a_celltype_key', n_meta_neighbors = 3)
  ```
  With `n_meta_neighbors` of 1, each cell type only has one neighboring cell type, that is, itself. This will result in strongly separated cell types in the final UMAP. Increasing `n_meta_neighbors` will loosen this restriction. For example, a `n_meta_neighbors` of 2 allows each cell type to have, in addition to itself, one nearest neighboring cell type based on the transcriptomic distances calculated by CellTypist. This parameter defaults to 3, meaning that a linear spectrum of transcriptomic structure can possibly exist for each cell type.
  </details>
</details>

<details>
<summary><strong>2. Tips for data integration</strong></summary>

+ <details>
  <summary><strong>2.1. Partial annotation</strong></summary>

  Partial annotation (an `.obs` column combining annotated and unannotated cells) is allowed as the `cell_type` parameter of `celltypist.integrate`. You need to explicitly name unannotated cells as `'UNASSIGNED'` for use in CellTypist (definition of symbols can be found [here](https://github.com/Teichlab/celltypist/blob/main/celltypist/contro/symbols.py)).  
    
  If low-hierarchy cell reannotation is used as input, it naturally contains `'UNASSIGNED'` cells (see the harmonisation section `2.2.`). These unannotated cells will be added to the search space of other cell types, and meanwhile will search neighbors against all cell types, resulting in less supervised integration and a prolonged runtime. Thus it is recommended to start with high-hierarchy cell types where all cells are fully annotated if you would like to feed the upstream harmonised cell types as input for the downstream integration.
  </details>

+ <details>
  <summary><strong>2.2. Rare cell types</strong></summary>

  When an abundant cell type is annotated/distributed across multiple batches (e.g., datasets), sometimes not all batches can harbour adequate numbers. This leads to a rare cell type defined within the context of a specific batch. During neighborhood construction, if this batch cannot provide enough neighboring cells for this cell type, search space will be expanded to all cells in this batch.  
    
  Although this represents a safe solution in CellTypist to anchor nearest neighbors for rare cell types, runtime of the algorithm will be increased and cells from this cell type may not be robustly clustered. Keeping them is fine for CellTypist, but you can also remove such rare cell types in associated batches before running `celltypist.integrate` (a cell type with only a small number in a given batch naturally means that this batch may not be qualified for hosting this cell type). Example code is:
  ```python
  #Remove cells from cell types that have <=5 cells in a batch.
  combined = adata.obs['a_batch_key'].astype(str) + adata.obs['a_celltype_key'].astype(str)
  combined_counts = combined.value_counts()
  remove_combn = combined_counts.index[combined_counts <= 5]
  adata = adata[~combined.isin(remove_combn)].copy()
  ```
  </details>

+ <details>
  <summary><strong>2.3. Use CellTypist models for annotation and integration</strong></summary>

  `celltypist.integrate` requires cell annotation to be stored in the AnnData. This information can be obtained by different means. One quick way is to use available CellTypist models to annotate the data of interest (see the CellTypist model list [here](https://www.celltypist.org/models)).
  ```python
  #Annotate the data with a relevant model (immune model as an example here).
  adata = celltypist.annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = True).to_adata()
  ```
  Then integrate cells on the basis of the predicted cell types.
  ```python
  #`cell_type` can also be 'majority_voting'.
  celltypist.integrate(adata, batch = 'a_batch_key', cell_type = 'predicted_labels')
  ```
  Even the model does not exactly match the data (e.g., using an immune model to annotate a lung data), this approach can be still useful as cells from the same cell type will probably be assigned the same identity by the model, therefore containing information with respect to which cells should be placed together in the neighborhood graph.
  </details>
</details>

# Citation
Xu et al., Automatic cell type harmonization and integration across Human Cell Atlas datasets. bioRxiv (2023). [Preprint](https://doi.org/10.1101/2023.05.01.538994)  
  
Dominguez Conde et al., Cross-tissue immune cell analysis reveals tissue-specific features in humans. Science 376, eabl5197 (2022). [Link](https://doi.org/10.1126/science.abl5197)
