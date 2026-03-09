import numpy as np
import pandas as pd
import scanpy as sc
import os
from anndata import AnnData
from sklearn.preprocessing import StandardScaler,LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import SGDClassifier
from sklearn import __version__ as skv
from typing import Optional, Union
from .models import Model, HierModel
from . import logger
from .tree import Tree
from scipy.sparse import spmatrix
from datetime import datetime
import sys
import copy
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

def _prepare_data(X, labels, genes, transpose, check_expression, indent) -> tuple:
    """
    For internal use. Prepare data for celltypist training.
    """
    if (X is None) or (labels is None):
        raise ValueError(
                "🛑 Missing training data and/or training labels. Please provide both arguments")
    if isinstance(X, AnnData) or (isinstance(X, str) and X.endswith('.h5ad')):
        adata = sc.read(X) if isinstance(X, str) else X
        adata.var_names_make_unique()
        if adata.X[:1000].min() < 0:
            logger.info(f"{indent}👀 Detected scaled expression in the input data, will try the .raw attribute")
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
                raise ValueError(
                        "🛑 Missing `genes`. Please provide this argument together with the input mtx file")
            genes = _to_vector(genes)
            if len(genes) != adata.n_vars:
                raise ValueError(
                        f"🛑 The number of genes provided does not match the number of genes in {X}")
            adata.var_names = np.asarray(genes)
        adata.var_names_make_unique()
        if not float(adata.X[:1000].max()).is_integer():
            logger.warn(f"{indent}⚠️ Warning: the input file seems not a raw count matrix. The trained model may be biased")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        indata = adata.X
        genes = adata.var_names
        labels = _to_vector(labels)
    elif isinstance(X, str):
        raise ValueError(
                "🛑 Invalid input. Supported types: .csv, .txt, .tsv, .tab, .mtx, .mtx.gz and .h5ad")
    else:
        #logger.info(f"{indent}👀 The input training data is processed as an array-like object")
        indata = X
        if transpose:
            indata = indata.transpose()
        if isinstance(indata, pd.DataFrame):
            genes = indata.columns
        else:
            if genes is None:
                raise ValueError(
                        "🛑 Missing `genes`. Please provide this argument together with the input training data")
            genes = _to_vector(genes)
        labels = _to_vector(labels)
    indata = indata if isinstance(indata, spmatrix) else np.asarray(indata)
    labels = np.asarray(labels)
    genes = np.asarray(genes)
    if check_expression and (np.abs(np.expm1(indata[0]).sum()-10000) > 1):
        raise ValueError(
                "🛑 Invalid expression matrix, expect log1p normalized expression to 10000 counts per cell")
    if len(labels) != indata.shape[0]:
        raise ValueError(
                f"🛑 Length of training labels ({len(labels)}) does not match the number of input cells ({indata.shape[0]})")
    if len(genes) != indata.shape[1]:
        raise ValueError(
                f"🛑 The number of genes ({len(genes)}) provided does not match the number of genes in the training data ({indata.shape[1]})")
    return indata, labels, genes

def _LRClassifier(indata, labels, C, solver, max_iter, n_jobs, indent, **kwargs) -> LogisticRegression:
    """
    For internal use. Get the logistic Classifier.
    """
    no_cells = len(labels)
    if solver is None:
        solver = 'sag' if no_cells>50000 else 'lbfgs'
    elif solver not in ('liblinear', 'lbfgs', 'newton-cg', 'sag', 'saga'):
        raise ValueError(
                f"🛑 Invalid `solver`, should be one of `'liblinear'`, `'lbfgs'`, `'newton-cg'`, `'sag'`, and `'saga'`")
    logger.info(f"{indent}🏋️ Training data using logistic regression")
    if (no_cells > 100000) and (indata.shape[1] > 10000):
        logger.warn(f"{indent}⚠️ Warning: it may take a long time to train this dataset with {no_cells} cells and {indata.shape[1]} genes, try to downsample cells and/or restrict genes to a subset (e.g., hvgs)")
    classifier = LogisticRegression(C = C, solver = solver, max_iter = max_iter, multi_class = 'ovr', n_jobs = n_jobs, **kwargs)
    classifier.fit(indata, labels)
    return classifier

def _cuLRClassifier(indata, labels, C, solver, max_iter, indent, **kwargs) -> LogisticRegression:
    """
    For internal use. Get the cuml logistic Classifier.
    """
    solver = 'qn' if solver is None else solver
    if solver != 'qn':
        raise ValueError(
                f"🛑 Invalid `solver`, should be `'qn'` to run on GPU")
    le = LabelEncoder()
    labels_ = le.fit_transform(labels)
    logger.info(f"{indent}🏋️ Training data using logistic regression on GPU")
    no_cells = len(labels)
    if (no_cells > 100000) and (indata.shape[1] > 10000):
        logger.warn(f"{indent}⚠️ Warning: it may take a long time to train this dataset with {no_cells} cells and {indata.shape[1]} genes, try to downsample cells and/or restrict genes to a subset (e.g., hvgs)")
    classifier_ = cuLogisticRegression(C = C, max_iter = max_iter, solver = solver, **kwargs)
    classifier_.fit(indata, labels_)
    classifier = LogisticRegression(multi_class = 'ovr')
    for attr in ('C', 'class_weight', 'fit_intercept', 'l1_ratio', 'max_iter', 'penalty', 'tol', 'solver', 'coef_', 'intercept_', 'verbose'):
        setattr(classifier, attr, getattr(classifier_, attr))
    classifier.classes_ = le.inverse_transform(classifier_.classes_)
    return classifier

def _SGDClassifier(indata, labels,
                   alpha, max_iter, n_jobs,
                   mini_batch, batch_number, batch_size, epochs, balance_cell_type, indent, **kwargs) -> SGDClassifier:
    """
    For internal use. Get the SGDClassifier.
    """
    loss_mode = 'log_loss' if float(skv[:3]) >= 1.1 else 'log'
    classifier = SGDClassifier(loss = loss_mode, alpha = alpha, max_iter = max_iter, n_jobs = n_jobs, **kwargs)
    if not mini_batch:
        logger.info(f"{indent}🏋️ Training data using SGD logistic regression")
        if (len(labels) > 100000) and (indata.shape[1] > 10000):
            logger.warn(f"{indent}⚠️ Warning: it may take a long time to train this dataset with {len(labels)} cells and {indata.shape[1]} genes, try to downsample cells and/or restrict genes to a subset (e.g., hvgs)")
        classifier.fit(indata, labels)
    else:
        logger.info(f"{indent}🏋️ Training data using mini-batch SGD logistic regression")
        no_cells = len(labels)
        if no_cells < 10000:
            logger.warn(f"{indent}⚠️ Warning: the number of cells ({no_cells}) is not big enough to conduct a proper mini-batch training. You may consider using traditional SGD classifier (mini_batch = False)")
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
        unique_labels = np.unique(labels)
        for epoch in range(1, (epochs+1)):
            logger.info(f"{indent}⏳ Epochs: [{epoch}/{epochs}]")
            if not balance_cell_type:
                sampled_cell_index = np.random.choice(no_cells, no_cells_sample, replace = False)
            else:
                sampled_cell_index = np.random.choice(no_cells, no_cells_sample, replace = False, p = p)
            for start in starts:
                s_index = sampled_cell_index[start:start+batch_size]
                classifier.partial_fit(indata[s_index], labels[s_index], classes = unique_labels)
    return classifier

def _prepare_params(X, labels, genes, transpose_input, with_mean, check_expression, max_iter, indent, copy) -> tuple:
    """
    For internal use. Wrapper code before the actual classifier.
    """
    #prepare
    logger.info(f"{indent}🍳 Preparing data before training")
    indata, labels, genes = _prepare_data(X, labels, genes, transpose_input, check_expression, indent)
    #filter
    if isinstance(indata, spmatrix):
        flag = indata.getnnz(axis = 0) == 0
    else:
        flag = np.count_nonzero(indata, axis = 0) == 0
    if flag.any():
        logger.info(f"{indent}✂️ {flag.sum()} non-expressed genes are filtered out")
        indata = indata[:, ~flag]
        genes = genes[~flag]
        copy = False
    #report data stats
    logger.info(f"{indent}🔬 Input data has {indata.shape[0]} cells and {indata.shape[1]} genes")
    if with_mean and isinstance(indata, spmatrix):
        indata = indata.toarray()
        copy = False
    #scaler
    logger.info(f"{indent}⚖️ Scaling input data")
    scaler = StandardScaler(with_mean = with_mean, copy = copy)
    indata = scaler.fit_transform(indata)
    if isinstance(indata, spmatrix):
        np.minimum(indata.data, 10, out = indata.data)
    else:
        np.minimum(indata, 10, out = indata)
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
    return indata, labels, genes, max_iter, scaler

def _actual_classifier(indata, labels, genes, max_iter, scaler,
        C, solver, n_jobs, use_SGD, alpha, use_GPU, mini_batch, batch_number, batch_size, epochs, balance_cell_type, feature_selection, top_genes, date, details, url, source, version, indent, **kwargs) -> Model:
    """
    For internal use. The actual classifier.
    """
    #classifier
    if use_SGD or feature_selection:
        classifier = _SGDClassifier(indata = indata, labels = labels, alpha = alpha, max_iter = max_iter, n_jobs = n_jobs, mini_batch = mini_batch, batch_number = batch_number, batch_size = batch_size, epochs = epochs, balance_cell_type = balance_cell_type, indent = indent, **kwargs)
    elif use_GPU:
        classifier = _cuLRClassifier(indata = indata, labels = labels, C = C, solver = solver, max_iter = max_iter, indent = indent, **kwargs)
    else:
        classifier = _LRClassifier(indata = indata, labels = labels, C = C, solver = solver, max_iter = max_iter, n_jobs = n_jobs, indent = indent, **kwargs)
    #feature selection -> new classifier and scaler
    if feature_selection:
        logger.info(f"{indent}🔎 Selecting features")
        if len(genes) <= top_genes:
            raise ValueError(
                    f"🛑 The number of genes ({len(genes)}) is fewer than the `top_genes` ({top_genes}). Unable to perform feature selection")
        gene_index = np.argpartition(np.abs(classifier.coef_), -top_genes, axis = 1)[:, -top_genes:]
        gene_index = np.unique(gene_index)
        logger.info(f"{indent}🧬 {len(gene_index)} features are selected")
        genes = genes[gene_index]
        indata = indata[:, gene_index]
        logger.info(f"{indent}🏋️ Starting the second round of training")
        if use_SGD:
            classifier = _SGDClassifier(indata = indata, labels = labels, alpha = alpha, max_iter = max_iter, n_jobs = n_jobs, mini_batch = mini_batch, batch_number = batch_number, batch_size = batch_size, epochs = epochs, balance_cell_type = balance_cell_type, indent = indent, **kwargs)
        elif use_GPU:
            classifier = _cuLRClassifier(indata = indata, labels = labels, C = C, solver = solver, max_iter = max_iter, indent = indent, **kwargs)
        else:
            classifier = _LRClassifier(indata = indata, labels = labels, C = C, solver = solver, max_iter = max_iter, n_jobs = n_jobs, indent = indent, **kwargs)
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
    logger.info(f"{indent}✅ Model training done!")
    return Model(classifier, scaler, description)

def train(X = None,
          labels: Optional[Union[str, list, tuple, np.ndarray, pd.Series, pd.Index]] = None,
          genes: Optional[Union[str, list, tuple, np.ndarray, pd.Series, pd.Index]] = None,
          transpose_input: bool = False,
          copy: bool = True, with_mean: bool = True,
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
    copy
        Whether to make a copy of input data for data scaling.
        (Default: `True`)
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
    #prepare params
    copy = False if isinstance(X, str) else copy
    indata, labels, genes, max_iter, scaler = _prepare_params(X, labels, genes, transpose_input, with_mean, check_expression, max_iter, '', copy)
    #actual classifier
    model = _actual_classifier(indata, labels, genes, max_iter, scaler, C, solver, n_jobs, use_SGD, alpha, use_GPU, mini_batch, batch_number, batch_size, epochs, balance_cell_type, feature_selection, top_genes, date, details, url, source, version, '', **kwargs)
    return model

def hier_train(X = None,
               tree: Optional[Union[Tree, str, dict]] = None,
               leaf_anno: Optional[Union[str, list, tuple, np.ndarray, pd.Series, pd.Index]] = None,
               genes: Optional[Union[str, list, tuple, np.ndarray, pd.Series, pd.Index]] = None,
               transpose_input: bool = False,
               copy: bool = True, with_mean: bool = True,
               check_expression: bool = True,
               mode: str = 'LCPN',
               C: float = 1.0, solver: Optional[str] = None, max_iter: Optional[int] = None, n_jobs: Optional[int] = None,
               use_SGD: bool = False, alpha: float = 0.0001,
               use_GPU: bool = False,
               mini_batch: bool = False, batch_number: int = 100, batch_size: int = 1000, epochs: int = 10, balance_cell_type: bool = False,
               feature_selection: bool = True, top_genes: int = 300,
               date: str = '', details: str = '', url: str = '', source: str = '', version: str = '',
               save_strategy: str = 'checkpointed', out_dir: Optional[str] = None, resume: bool = False,
               **kwargs) -> HierModel:
    """
    Train a hierarchical ensemble model for cell type classification, built from local classifiers using either LCPN (Local Classifier per Parent Node) or LCL (Local Classifier per Level).

    Parameters
    ----------
    X
        Path to the input count matrix (supported types are csv, txt, tsv, tab and mtx) or AnnData (h5ad).
        Also accepts the input as an :class:`~anndata.AnnData` object, or any array-like objects already loaded in memory.
        See `check_expression` for detailed format requirements.
        A cell-by-gene format is desirable (see `transpose_input` for more information).
    tree
        A :class:`~celltypist.tree.Tree` object representing the predefined cell type hierarchy.
        Also accepts the tree in dictionary form, or via a JSON file with that representation.
        This argument is ignored if a previous run is resumed (`save_strategy = 'checkpointed'` and `resume = True`).
    leaf_anno
        Path to the file containing leaf cell type label per line corresponding to the cells in `X`.
        Also accepts any list-like objects already loaded in memory (such as an array).
        If `X` is specified as an AnnData, this argument can also be set as a column name from cell metadata.
    genes
        Path to the file containing one gene per line corresponding to the genes in `X`.
        Also accepts any list-like objects already loaded in memory (such as an array).
        Note `genes` will be extracted from `X` where possible (e.g., `X` is an AnnData or data frame).
    transpose_input
        Whether to transpose the input matrix. Set to `True` if `X` is provided in a gene-by-cell format.
        (Default: `False`)
    copy
        Whether to make a copy of input data for data scaling.
        (Default: `True`)
    with_mean
        Whether to subtract the mean values during data scaling. Setting to `False` can lower the memory usage when the input is a sparse matrix but may slightly reduce the model performance.
        (Default: `True`)
    check_expression
        Check whether the expression matrix in the input data is supplied as required.
        Except the case where a path to the raw count table file is specified, all other inputs for `X` should be in log1p normalized expression to 10000 counts per cell.
        Set to `False` if you want to train the data regardless of the expression formats.
        (Default: `True`)
    mode
        Local classifier to use, either `'LCPN'` or `'LCL'`.
        (Default: `'LCPN'`)
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
        (Default: `True`)
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
    save_strategy
        Runs to completion (`save_strategy = 'atomic'`) or saves intermediate files incrementally to `out_dir` (`save_strategy = 'checkpointed'`), before returning a :class:`~celltypist.models.HierModel` instance.
        (Default: `'checkpointed'`)
    out_dir
        Path to the model’s working directory where all local classifiers and the tree are written.
        This argument is only relevant if saving the model on-the-fly (`save_strategy = 'checkpointed'`).
        Default to `{tree.handle}_{mode}` in the current working directory if not provided.
    resume
        Whether to resume from an existing run in `out_dir`.
        This argument is only relevant if saving the model on-the-fly (`save_strategy = 'checkpointed'`).
        (Default: `False`)
    **kwargs
        Other keyword arguments passed to :class:`~sklearn.linear_model.LogisticRegression` (`use_SGD = False` and `use_GPU = False`), :class:`cuml.LogisticRegression` (`use_SGD = False` and `use_GPU = True`), or :class:`~sklearn.linear_model.SGDClassifier` (`use_SGD = True`).

    Returns
    ----------
    :class:`~celltypist.models.HierModel`
        A :class:`~celltypist.models.HierModel` object trained by celltypist.
    """
    #validate params
    copy = False if isinstance(X, str) else copy
    if not use_SGD and use_GPU and 'cuml' not in sys.modules:
        logger.warn(f"⚠️ Warning: to run logistic regression on GPU, please first install cuml")
        return
    if mode not in ('LCPN', 'LCL'):
        raise ValueError(
                f"🛑 Unrecognized `mode` value, should be one of `'LCPN'` or `'LCL'`")
    if save_strategy not in ('checkpointed', 'atomic'):
        raise ValueError(
                f"🛑 Unrecognized `save_strategy` value, should be one of `'checkpointed'` or `'atomic'`")
    #checkpointed: continued=True & loaded tree & out_dir (resume=True) vs. continued=False & empty out_dir (resume=False)
    #atomic: only continued=False
    continued = False
    if save_strategy == 'checkpointed':
        if out_dir is None:
            if resume:
                raise ValueError(
                        f"🛑 Please provide `out_dir` to resume a previous training run")
            else:
                out_dir = f"{tree.handle if isinstance(tree, Tree) else Tree.from_json(tree).handle}_{mode}"
                logger.info(f"📂 Output directory not specified. Using default: `{out_dir}`")
        if resume:
            tree_file = os.path.join(out_dir, "tree.json")
            if os.path.isfile(tree_file):
                tree = Tree.from_json(tree_file)
                if not hasattr(tree, 'mode'):
                    raise ValueError(
                            f"🛑 Cannot resume training: `{out_dir}` does not contain a previous run")
                if tree.mode != mode:
                    raise ValueError(
                            f"🛑 The mode of previous run in `{out_dir}` ('{tree.mode}') does not match the current training mode ('{mode}'). Please ensure you are resuming with the same mode")
                logger.info(f"📂 Resuming previous training run in `{out_dir}`")
                continued = True
            else:
                raise ValueError(
                        f"🛑 Cannot resume training. Make sure `out_dir` exists and contains a previous run")
        else:
            if os.path.isdir(out_dir):
                if os.listdir(out_dir):
                    raise ValueError(
                            f"🛑 Invalid output directory `{out_dir}`. Please remove its contents or specify an empty folder to start a new training run")
                else:
                    logger.info(f"📂 Starting a new training run in `{out_dir}`")
            else:
                os.mkdir(out_dir)
                logger.info(f"📂 Created new output directory `{out_dir}`. Starting a new training run")
    else:
        logger.info("⚛️ Using atomic save strategy; training will run to completion")
    #now all have continued, loaded tree, and out_dir (if needed) & model_mapping (even empty)
    model_mapping = {}
    if not continued:
        if tree is None:
            raise ValueError(
                    f"🛑 Please provide the `tree` argument")
        tree = tree.copy() if isinstance(tree, Tree) else Tree.from_json(tree)
        for node in tree.iter_nodes(leaf_only = False):
            node.model = ''
        attr_list = list(tree.__dict__)
        for attr in attr_list:
            if attr.startswith('level') and attr.endswith('_classifier'):
                delattr(tree, attr)
    else:
        logger.info(f"⏳ Loading previous run")
        if tree.mode == "LCPN":
            for node in tree.iter_nodes(leaf_only = False):
                if node.model:
                    model_mapping[node.model] = Model.load(os.path.join(out_dir, node.model))
        else:
            for attr, val in tree.__dict__.items():
                if attr.startswith("level") and attr.endswith("_classifier"):
                    model_mapping[val] = Model.load(os.path.join(out_dir, val))
    depth = tree.depth
    #real leaf_anno -> multi_anno & set size
    if isinstance(X, AnnData) or (isinstance(X, str) and X.endswith('.h5ad')):
        adata = sc.read(X, backed = 'r') if isinstance(X, str) else X
        if isinstance(leaf_anno, str) and (leaf_anno in adata.obs):
            leaf_anno = adata.obs[leaf_anno]
        else:
            leaf_anno = _to_vector(leaf_anno)
    else:
        leaf_anno = _to_vector(leaf_anno)
    leaf_anno = np.asarray(leaf_anno)
    multi_anno = tree.get_multilevel_anno(leaf_anno)
    if not continued:
        tree.assign_size(leaf_anno)
    else:
        proposed = tree.copy()
        proposed.assign_size(leaf_anno)
        p_sizes = [node.size for node in proposed.iter_nodes(leaf_only = False)]
        t_sizes = [node.size for node in tree.iter_nodes(leaf_only = False)]
        if not np.array_equal(p_sizes, t_sizes):
            raise ValueError(
                    f"🛑 The current `leaf_anno` does not match the one used in the previous run. Please resume with the same `leaf_anno`")
    #early return
    if mode == "LCPN":
        n_needed_models = 0
        for node in tree.iter_nodes(leaf_only = False):
            if sum(child.size > 0 for child in node.children) >= 2:
                n_needed_models += 1
    else:
        n_needed_models = (multi_anno.nunique(axis = 0) >= 2).sum()
    if continued and len(model_mapping) == n_needed_models:
        logger.info(f"✅ No need to resume, training in `{out_dir}` is already complete. The model is now loaded")
        return HierModel(tree, model_mapping, mode = mode, date = tree.date)
    #main
    if mode == 'LCL':
        indata, _, genes, max_iter, scaler = _prepare_params(X, leaf_anno, genes, transpose_input, with_mean, check_expression, max_iter, '', copy)
        logger.info(f"📚 Total models to train: {n_needed_models}")
        #LCL
        ith = 0
        for n in range(2, depth+1):
            labels = np.asarray(multi_anno[f"level{n}_anno"])
            if len(np.unique(labels)) < 2:
                continue
            ith += 1
            filename = f"{tree.handle}_level{n}.pkl"
            if filename in model_mapping:
                logger.info(f"⏩ Skipping level-{n} model training [{ith}/{n_needed_models}]: `{filename}` (model exists)")
                continue
            logger.info(f"🏋️ Training level-{n} model [{ith}/{n_needed_models}]: `{filename}`")
            model = _actual_classifier(indata, labels, genes, max_iter, copy.deepcopy(scaler), C, solver, n_jobs, use_SGD, alpha, use_GPU, mini_batch, batch_number, batch_size, epochs, balance_cell_type, feature_selection, top_genes, date, f"{details} (level {n})" if details else '', 'N/A', source, version, '      ', **kwargs)
            setattr(tree, f"level{n}_classifier", filename)
            model_mapping[filename] = model
            if save_strategy == 'checkpointed':
                hier_model = HierModel(tree, model_mapping, mode = mode, date = date, details = details, url = url, source = source, version = version)
                model.write(os.path.join(out_dir, filename))
                hier_model.tree.write(os.path.join(out_dir, 'tree.json'))
    else:
        logger.info(f"📚 Total models to train: {n_needed_models}")
        #Get subsettable X
        if isinstance(X, AnnData) or (isinstance(X, str) and X.endswith('.h5ad')):
            X = sc.read(X) if isinstance(X, str) else X
            X.var_names_make_unique()
            if X.X[:1000].min() < 0:
                logger.info(f"👀 Detected scaled expression in the input data, will try the .raw attribute")
                try:
                    genes = X.raw.var_names
                    X = X.raw.X
                except Exception as e:
                    raise Exception(
                            f"🛑 Fail to use the .raw attribute in the input object. {e}")
            else:
                genes = X.var_names
                X = X.X
            transpose_input = False
        elif isinstance(X, str) and X.endswith(('.csv', '.txt', '.tsv', '.tab', '.mtx', '.mtx.gz')):
            X_old = X
            X = sc.read(X)
            if transpose_input:
                X = X.transpose()
            if X_old.endswith(('.mtx', '.mtx.gz')):
                if genes is None:
                    raise ValueError(
                            "🛑 Missing `genes`. Please provide this argument together with the input mtx file")
                genes = _to_vector(genes)
                if len(genes) != X.n_vars:
                    raise ValueError(
                            f"🛑 The number of genes provided does not match the number of genes in {X_old}")
                X.var_names = np.asarray(genes)
            if not float(X.X[:1000].max()).is_integer():
                logger.warn(f"⚠️ Warning: the input file seems not a raw count matrix. The trained model may be biased")
            sc.pp.normalize_total(X, target_sum = 1e4)
            sc.pp.log1p(X)
            genes = X.var_names
            X = X.X
            transpose_input = False
        elif isinstance(X, str):
            raise ValueError(
                    "🛑 Invalid input. Supported types: .csv, .txt, .tsv, .tab, .mtx, .mtx.gz and .h5ad")
        else:
            if transpose_input:
                X = X.transpose()
                transpose_input = False
        #LCPN
        ith = 0
        for node in tree.iter_nodes(leaf_only = False):
            if sum(child.size > 0 for child in node.children) < 2:
                continue
            ith += 1
            filename = f"{node.internal_name}.pkl"
            if filename in model_mapping:
                logger.info(f"⏩ Skipping model training for node '{node.original_name}' [{ith}/{n_needed_models}]: `{filename}` (model exists)")
                continue
            node_depth = len(tree.extract_path(node.original_name, print_path = False))
            flag = (multi_anno[f"level{node_depth}_anno"] == node.original_name).values
            logger.info(f"🏋️ Training local model for node '{node.original_name}' [{ith}/{n_needed_models}]: `{filename}`")
            indata, labels, out_genes, out_max_iter, scaler = _prepare_params(X[flag], multi_anno[f"level{node_depth+1}_anno"][flag], genes, transpose_input, with_mean, check_expression, max_iter, '      ', copy)
            model = _actual_classifier(indata, labels, out_genes, out_max_iter, scaler, C, solver, n_jobs, use_SGD, alpha, use_GPU, mini_batch, batch_number, batch_size, epochs, balance_cell_type, feature_selection, top_genes, date, f"cell subtypes of {node.original_name}", 'N/A', source, version, '      ', **kwargs)
            setattr(node, 'model', filename)
            model_mapping[filename] = model
            if save_strategy == 'checkpointed':
                hier_model = HierModel(tree, model_mapping, mode = mode, date = date, details = details, url = url, source = source, version = version)
                model.write(os.path.join(out_dir, filename))
                hier_model.tree.write(os.path.join(out_dir, 'tree.json'))
    #done
    logger.info(f"✅ Hierarchical training completed successfully (mode = {mode})")
    return HierModel(tree, model_mapping, mode = mode, date = date, details = details, url = url, source = source, version = version)
