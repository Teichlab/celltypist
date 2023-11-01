import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from sklearn.preprocessing import StandardScaler,LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import SGDClassifier
from sklearn import __version__ as skv
from typing import Optional, Union
from .models import Model
from . import logger
from scipy.sparse import spmatrix
from datetime import datetime
import sys
try:
    from cuml import LogisticRegression as cuLogisticRegression
except ImportError:
    pass

def _to_vector(_vector_or_file):
    """
    For internal use. Turn a file into an array.
    """
    if isinstance(_vector_or_file, str):
        try:
            return pd.read_csv(_vector_or_file, header=None)[0].values
        except Exception as e:
            raise Exception(
                    f"🛑 {e}")
    else:
        return _vector_or_file

def _to_array(_array_like) -> np.ndarray:
    """
    For internal use. Turn an array-like object into an array.
    """
    if isinstance(_array_like, pd.DataFrame):
        return _array_like.values
    elif isinstance(_array_like, spmatrix):
        return _array_like.toarray()
    elif isinstance(_array_like, np.matrix):
        return np.array(_array_like)
    elif isinstance(_array_like, np.ndarray):
        return _array_like
    else:
        raise TypeError(
                f"🛑 Please provide a valid array-like object as input")

def _prepare_data(X, labels, genes, transpose) -> tuple:
    """
    For internal use. Prepare data for celltypist training.
    """
    if (X is None) or (labels is None):
        raise Exception(
                "🛑 Missing training data and/or training labels. Please provide both arguments")
    if isinstance(X, AnnData) or (isinstance(X, str) and X.endswith('.h5ad')):
        adata = sc.read(X) if isinstance(X, str) else X
        adata.var_names_make_unique()
        if adata.X[:1000].min() < 0:
            logger.info("👀 Detected scaled expression in the input data, will try the .raw attribute")
            try:
                indata = adata.raw.X
                genes = adata.raw.var_names
            except Exception as e:
                raise Exception(
                        f"🛑 Fail to use the .raw attribute in the input object. {e}")
        else:
            indata = adata.X
            genes = adata.var_names
        if isinstance(labels, str) and (labels in adata.obs):
            labels = adata.obs[labels]
        else:
            labels = _to_vector(labels)
    elif isinstance(X, str) and X.endswith(('.csv', '.txt', '.tsv', '.tab', '.mtx', '.mtx.gz')):
        adata = sc.read(X)
        if transpose:
            adata = adata.transpose()
        if X.endswith(('.mtx', '.mtx.gz')):
            if genes is None:
                raise Exception(
                        "🛑 Missing `genes`. Please provide this argument together with the input mtx file")
            genes = _to_vector(genes)
            if len(genes) != adata.n_vars:
                raise ValueError(
                        f"🛑 The number of genes provided does not match the number of genes in {X}")
            adata.var_names = np.array(genes)
        adata.var_names_make_unique()
        if not float(adata.X[:1000].max()).is_integer():
            logger.warn(f"⚠️ Warning: the input file seems not a raw count matrix. The trained model may be biased")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        indata = adata.X
        genes = adata.var_names
        labels = _to_vector(labels)
    elif isinstance(X, str):
        raise ValueError(
                "🛑 Invalid input. Supported types: .csv, .txt, .tsv, .tab, .mtx, .mtx.gz and .h5ad")
    else:
        logger.info("👀 The input training data is processed as an array-like object")
        indata = X
        if transpose:
            indata = indata.transpose()
        if isinstance(indata, pd.DataFrame):
            genes = indata.columns
        else:
            if genes is None:
                raise Exception(
                        "🛑 Missing `genes`. Please provide this argument together with the input training data")
            genes = _to_vector(genes)
        labels = _to_vector(labels)
    return indata, labels, genes

def _LRClassifier(indata, labels, C, solver, max_iter, n_jobs, **kwargs) -> LogisticRegression:
    """
    For internal use. Get the logistic Classifier.
    """
    no_cells = len(labels)
    if solver is None:
        solver = 'sag' if no_cells>50000 else 'lbfgs'
    elif solver not in ('liblinear', 'lbfgs', 'newton-cg', 'sag', 'saga'):
        raise ValueError(
                f"🛑 Invalid `solver`, should be one of `'liblinear'`, `'lbfgs'`, `'newton-cg'`, `'sag'`, and `'saga'`")
    logger.info(f"🏋️ Training data using logistic regression")
    if (no_cells > 100000) and (indata.shape[1] > 10000):
        logger.warn(f"⚠️ Warning: it may take a long time to train this dataset with {no_cells} cells and {indata.shape[1]} genes, try to downsample cells and/or restrict genes to a subset (e.g., hvgs)")
    classifier = LogisticRegression(C = C, solver = solver, max_iter = max_iter, multi_class = 'ovr', n_jobs = n_jobs, **kwargs)
    classifier.fit(indata, labels)
    return classifier

def _cuLRClassifier(indata, labels, C, solver, max_iter, **kwargs) -> LogisticRegression:
    """
    For internal use. Get the cuml logistic Classifier.
    """
    solver = 'qn' if solver is None else solver
    if solver != 'qn':
        raise ValueError(
                f"🛑 Invalid `solver`, should be `'qn'` to run on GPU")
    le = LabelEncoder()
    labels_ = le.fit_transform(labels)
    logger.info(f"🏋️ Training data using logistic regression on GPU")
    no_cells = len(labels)
    if (no_cells > 100000) and (indata.shape[1] > 10000):
        logger.warn(f"⚠️ Warning: it may take a long time to train this dataset with {no_cells} cells and {indata.shape[1]} genes, try to downsample cells and/or restrict genes to a subset (e.g., hvgs)")
    classifier_ = cuLogisticRegression(C = C, max_iter = max_iter, solver = solver, **kwargs)
    classifier_.fit(indata, labels_)
    classifier = LogisticRegression(multi_class = 'ovr')
    for attr in ['C', 'class_weight', 'fit_intercept', 'l1_ratio', 'max_iter', 'penalty', 'tol', 'solver', 'coef_', 'intercept_', 'verbose']:
        setattr(classifier, attr, getattr(classifier_, attr))
    classifier.classes_ = le.inverse_transform(classifier_.classes_)
    return classifier

def _SGDClassifier(indata, labels,
                   alpha, max_iter, n_jobs,
                   mini_batch, batch_number, batch_size, epochs, balance_cell_type, **kwargs) -> SGDClassifier:
    """
    For internal use. Get the SGDClassifier.
    """
    loss_mode = 'log_loss' if float(skv[:3]) >= 1.1 else 'log'
    classifier = SGDClassifier(loss = loss_mode, alpha = alpha, max_iter = max_iter, n_jobs = n_jobs, **kwargs)
    if not mini_batch:
        logger.info(f"🏋️ Training data using SGD logistic regression")
        if (len(labels) > 100000) and (indata.shape[1] > 10000):
            logger.warn(f"⚠️ Warning: it may take a long time to train this dataset with {len(labels)} cells and {indata.shape[1]} genes, try to downsample cells and/or restrict genes to a subset (e.g., hvgs)")
        classifier.fit(indata, labels)
    else:
        logger.info(f"🏋️ Training data using mini-batch SGD logistic regression")
        no_cells = len(labels)
        if no_cells < 10000:
            logger.warn(f"⚠️ Warning: the number of cells ({no_cells}) is not big enough to conduct a proper mini-batch training. You may consider using traditional SGD classifier (mini_batch = False)")
        if no_cells <= batch_size:
            raise ValueError(
                    f"🛑 Number of cells ({no_cells}) is fewer than the batch size ({batch_size}). Decrease `batch_size`, or use SGD directly (mini_batch = False)")
        no_cells_sample = min([batch_number*batch_size, no_cells])
        starts = np.arange(0, no_cells_sample, batch_size)
        if balance_cell_type:
            celltype_freq = np.unique(labels, return_counts = True)
            len_celltype = len(celltype_freq[0])
            mapping = pd.Series(1 / (celltype_freq[1]*len_celltype), index = celltype_freq[0])
            p = mapping[labels].values
        for epoch in range(1, (epochs+1)):
            logger.info(f"⏳ Epochs: [{epoch}/{epochs}]")
            if not balance_cell_type:
                sampled_cell_index = np.random.choice(no_cells, no_cells_sample, replace = False)
            else:
                sampled_cell_index = np.random.choice(no_cells, no_cells_sample, replace = False, p = p)
            for start in starts:
                classifier.partial_fit(indata[sampled_cell_index[start:start+batch_size]], labels[sampled_cell_index[start:start+batch_size]], classes = np.unique(labels))
    return classifier

def train(X = None,
          labels: Optional[Union[str, list, tuple, np.ndarray, pd.Series, pd.Index]] = None,
          genes: Optional[Union[str, list, tuple, np.ndarray, pd.Series, pd.Index]] = None,
          transpose_input: bool = False,
          with_mean: bool = True,
          check_expression: bool = True,
          #LR param
          C: float = 1.0, solver: Optional[str] = None, max_iter: Optional[int] = None, n_jobs: Optional[int] = None,
          #SGD param
          use_SGD: bool = False, alpha: float = 0.0001,
          #GPU param
          use_GPU: bool = False,
          #mini-batch
          mini_batch: bool = False, batch_number: int = 100, batch_size: int = 1000, epochs: int = 10, balance_cell_type: bool = False,
          #feature selection
          feature_selection: bool = False, top_genes: int = 300,
          #description
          date: str = '', details: str = '', url: str = '', source: str = '', version: str = '',
          #other param
          **kwargs
         ) -> Model:
    """
    Train a celltypist model using mini-batch (optional) logistic classifier with a global solver or stochastic gradient descent (SGD) learning.

    Parameters
    ----------
    X
        Path to the input count matrix (supported types are csv, txt, tsv, tab and mtx) or AnnData (h5ad).
        Also accepts the input as an :class:`~anndata.AnnData` object, or any array-like objects already loaded in memory.
        See `check_expression` for detailed format requirements.
        A cell-by-gene format is desirable (see `transpose_input` for more information).
    labels
        Path to the file containing cell type label per line corresponding to the cells in `X`.
        Also accepts any list-like objects already loaded in memory (such as an array).
        If `X` is specified as an AnnData, this argument can also be set as a column name from cell metadata.
    genes
        Path to the file containing one gene per line corresponding to the genes in `X`.
        Also accepts any list-like objects already loaded in memory (such as an array).
        Note `genes` will be extracted from `X` where possible (e.g., `X` is an AnnData or data frame).
    transpose_input
        Whether to transpose the input matrix. Set to `True` if `X` is provided in a gene-by-cell format.
        (Default: `False`)
    with_mean
        Whether to subtract the mean values during data scaling. Setting to `False` can lower the memory usage when the input is a sparse matrix but may slightly reduce the model performance.
        (Default: `True`)
    check_expression
        Check whether the expression matrix in the input data is supplied as required.
        Except the case where a path to the raw count table file is specified, all other inputs for `X` should be in log1p normalized expression to 10000 counts per cell.
        Set to `False` if you want to train the data regardless of the expression formats.
        (Default: `True`)
    C
        Inverse of L2 regularization strength for traditional logistic classifier. A smaller value can possibly improve model generalization while at the cost of decreased accuracy.
        This argument is ignored if SGD learning is enabled (`use_SGD = True`).
        (Default: 1.0)
    solver
        Algorithm to use in the optimization problem for traditional logistic classifier.
        The default behavior is to choose the solver according to the size of the input data.
        This argument is ignored if SGD learning is enabled (`use_SGD = True`).
    max_iter
        Maximum number of iterations before reaching the minimum of the cost function.
        Try to decrease `max_iter` if the cost function does not converge for a long time.
        This argument is for both traditional and SGD logistic classifiers, and will be ignored if mini-batch SGD training is conducted (`use_SGD = True` and `mini_batch = True`).
        Default to 200, 500, and 1000 for large (>500k cells), medium (50-500k), and small (<50k) datasets, respectively.
    n_jobs
        Number of CPUs used. Default to one CPU. `-1` means all CPUs are used.
        This argument is for both traditional and SGD logistic classifiers.
    use_SGD
        Whether to implement SGD learning for the logistic classifier.
        (Default: `False`)
    alpha
        L2 regularization strength for SGD logistic classifier. A larger value can possibly improve model generalization while at the cost of decreased accuracy.
        This argument is ignored if SGD learning is disabled (`use_SGD = False`).
        (Default: 0.0001)
    use_GPU
        Whether to use GPU for logistic classifier.
        This argument is ignored if SGD learning is enabled (`use_SGD = True`).
        (Default: `False`)
    mini_batch
        Whether to implement mini-batch training for the SGD logistic classifier.
        Setting to `True` may improve the training efficiency for large datasets (for example, >100k cells).
        This argument is ignored if SGD learning is disabled (`use_SGD = False`).
        (Default: `False`)
    batch_number
        The number of batches used for training in each epoch. Each batch contains `batch_size` cells.
        For datasets which cannot be binned into `batch_number` batches, all batches will be used.
        This argument is relevant only if mini-batch SGD training is conducted (`use_SGD = True` and `mini_batch = True`).
        (Default: 100)
    batch_size
        The number of cells within each batch.
        This argument is relevant only if mini-batch SGD training is conducted (`use_SGD = True` and `mini_batch = True`).
        (Default: 1000)
    epochs
        The number of epochs for the mini-batch training procedure.
        The default values of `batch_number`, `batch_size`, and `epochs` together allow observing ~10^6 training cells.
        This argument is relevant only if mini-batch SGD training is conducted (`use_SGD = True` and `mini_batch = True`).
        (Default: 10)
    balance_cell_type
        Whether to balance the cell type frequencies in mini-batches during each epoch.
        Setting to `True` will sample rare cell types with a higher probability, ensuring close-to-even cell type distributions in mini-batches.
        This argument is relevant only if mini-batch SGD training is conducted (`use_SGD = True` and `mini_batch = True`).
        (Default: `False`)
    feature_selection
        Whether to perform two-pass data training where the first round is used for selecting important features/genes using SGD learning.
        If `True`, the training time will be longer.
        (Default: `False`)
    top_genes
        The number of top genes selected from each class/cell-type based on their absolute regression coefficients.
        The final feature set is combined across all classes (i.e., union).
        (Default: 300)
    date
        Free text of the date of the model. Default to the time when the training is completed.
    details
        Free text of the description of the model.
    url
        Free text of the (possible) download url of the model.
    source
        Free text of the source (publication, database, etc.) of the model.
    version
        Free text of the version of the model.
    **kwargs
        Other keyword arguments passed to :class:`~sklearn.linear_model.LogisticRegression` (`use_SGD = False` and `use_GPU = False`), :class:`cuml.LogisticRegression` (`use_SGD = False` and `use_GPU = True`), or :class:`~sklearn.linear_model.SGDClassifier` (`use_SGD = True`).

    Returns
    ----------
    :class:`~celltypist.models.Model`
        An instance of the :class:`~celltypist.models.Model` trained by celltypist.
    """
    #Test GPU
    if not use_SGD and use_GPU and 'cuml' not in sys.modules:
        logger.warn(f"⚠️ Warning: to run logistic regression on GPU, please first install cuml")
        return
    #prepare
    logger.info("🍳 Preparing data before training")
    indata, labels, genes = _prepare_data(X, labels, genes, transpose_input)
    if isinstance(indata, pd.DataFrame):
        indata = indata.values
    elif with_mean and isinstance(indata, spmatrix):
        indata = indata.toarray()
    labels = np.array(labels)
    genes = np.array(genes)
    #check
    if check_expression and (np.abs(np.expm1(indata[0]).sum()-10000) > 1):
        raise ValueError(
                "🛑 Invalid expression matrix, expect log1p normalized expression to 10000 counts per cell")
    if len(labels) != indata.shape[0]:
        raise ValueError(
                f"🛑 Length of training labels ({len(labels)}) does not match the number of input cells ({indata.shape[0]})")
    if len(genes) != indata.shape[1]:
        raise ValueError(
                f"🛑 The number of genes ({len(genes)}) provided does not match the number of genes in the training data ({indata.shape[1]})")
    #filter
    flag = indata.sum(axis = 0) == 0
    if isinstance(flag, np.matrix):
        flag = flag.A1
    if flag.sum() > 0:
        logger.info(f"✂️ {flag.sum()} non-expressed genes are filtered out")
        #indata = indata[:, ~flag]
        genes = genes[~flag]
    #report data stats
    logger.info(f"🔬 Input data has {indata.shape[0]} cells and {(~flag).sum()} genes")
    #scaler
    logger.info(f"⚖️ Scaling input data")
    scaler = StandardScaler(with_mean = with_mean)
    indata = scaler.fit_transform(indata[:, ~flag] if flag.sum() > 0 else indata)
    indata[indata > 10] = 10
    #sklearn (Cython) does not support very large sparse matrices for the time being
    if isinstance(indata, spmatrix) and ((indata.indices.dtype == 'int64') or (indata.indptr.dtype == 'int64')):
        indata = indata.toarray()
    #max_iter
    if max_iter is None:
        if indata.shape[0] < 50000:
            max_iter = 1000
        elif indata.shape[0] < 500000:
            max_iter = 500
        else:
            max_iter = 200
    #classifier
    if use_SGD or feature_selection:
        classifier = _SGDClassifier(indata = indata, labels = labels, alpha = alpha, max_iter = max_iter, n_jobs = n_jobs, mini_batch = mini_batch, batch_number = batch_number, batch_size = batch_size, epochs = epochs, balance_cell_type = balance_cell_type, **kwargs)
    elif use_GPU:
        classifier = _cuLRClassifier(indata = indata, labels = labels, C = C, solver = solver, max_iter = max_iter, **kwargs)
    else:
        classifier = _LRClassifier(indata = indata, labels = labels, C = C, solver = solver, max_iter = max_iter, n_jobs = n_jobs, **kwargs)
    #feature selection -> new classifier and scaler
    if feature_selection:
        logger.info(f"🔎 Selecting features")
        if len(genes) <= top_genes:
            raise ValueError(
                    f"🛑 The number of genes ({len(genes)}) is fewer than the `top_genes` ({top_genes}). Unable to perform feature selection")
        gene_index = np.argpartition(np.abs(classifier.coef_), -top_genes, axis = 1)[:, -top_genes:]
        gene_index = np.unique(gene_index)
        logger.info(f"🧬 {len(gene_index)} features are selected")
        genes = genes[gene_index]
        #indata = indata[:, gene_index]
        logger.info(f"🏋️ Starting the second round of training")
        if use_SGD:
            classifier = _SGDClassifier(indata = indata[:, gene_index], labels = labels, alpha = alpha, max_iter = max_iter, n_jobs = n_jobs, mini_batch = mini_batch, batch_number = batch_number, batch_size = batch_size, epochs = epochs, balance_cell_type = balance_cell_type, **kwargs)
        elif use_GPU:
            classifier = _cuLRClassifier(indata = indata[:, gene_index], labels = labels, C = C, solver = solver, max_iter = max_iter, **kwargs)
        else:
            classifier = _LRClassifier(indata = indata[:, gene_index], labels = labels, C = C, solver = solver, max_iter = max_iter, n_jobs = n_jobs, **kwargs)
        scaler.mean_ = scaler.mean_[gene_index]
        scaler.var_ = scaler.var_[gene_index]
        scaler.scale_ = scaler.scale_[gene_index]
        scaler.n_features_in_ = len(gene_index)
    #model finalization
    classifier.features = genes
    classifier.n_features_in_ = len(genes)
    if not date:
        date = str(datetime.now())
    description = {'date': date, 'details': details, 'url': url, 'source': source, 'version': version, 'number_celltypes': len(classifier.classes_)}
    logger.info(f"✅ Model training done!")
    return Model(classifier, scaler, description)
