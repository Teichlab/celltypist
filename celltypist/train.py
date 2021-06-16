import numpy as np
import pandas as pd
import scanpy as sc
from scanpy import AnnData
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import SGDClassifier
from typing import Optional, Union
from .models import Model
from . import logger
from scipy.sparse import spmatrix

def _load_file_as_list(_file):
    """
    For internal use. Load file if it is a string.
    """
    if isinstance(_file, str):
        try:
            with open(_file, 'rt') as f:
                return [x.strip() for x in f.readlines()]
        except Exception as e:
            raise Exception(f"🛑 {e}")
    else:
        return _file

def _prepare_data(X, labels, genes, transpose):
    """
    For internal use. Prepare data for celltypist training.
    """
    if (X is None) or (labels is None):
        raise Exception("🛑 Missing training data and/or training labels. Please provide both arguments")
    if isinstance(X, AnnData) or (isinstance(X, str) and X.endswith('.h5ad')):
        adata = sc.read(X) if isinstance(X, str) else X
        if adata.X.min() < 0:
            logger.info("👀 Detected scaled expression in the input data, will try the .raw attribute")
            try:
                indata = adata.raw.X.copy()
                genes = adata.raw.var_names.copy()
            except Exception as e:
                raise Exception(f"🛑 Fail to use the .raw attribute in the input object. {e}")
        else:
            indata = adata.X.copy()
            genes = adata.var_names.copy()
        if np.abs(np.expm1(indata[0]).sum()-10000) > 1:
            raise ValueError("🛑 Invalid expression matrix, expect log1p normalized expression to 10000 counts per cell")
        if isinstance(labels, str):
            if labels in adata.obs:
                labels = adata.obs[labels]
            else:
                labels = _load_file_as_list(labels)
        if len(labels) != indata.shape[0]:
            raise ValueError(f"🛑 Length of training labels ({len(labels)}) does not match the number of input cells ({indata.shape[0]})")
    elif isinstance(X, str) and X.endswith(('.csv', '.txt', '.tsv', '.tab', '.mtx', '.mtx.gz')):
        adata = sc.read(X)
        if transpose:
            adata = adata.transpose()
        if X.endswith(('.mtx', '.mtx.gz')):
            if genes is None:
                raise Exception("🛑 Missing `genes`. Please provide this argument together with the input mtx file")
            genes = _load_file_as_list(genes)
            if len(genes) != adata.n_vars:
                raise ValueError(f"🛑 The number of genes does not match the number of genes in {X}")
            adata.var_names = genes
        adata.var_names_make_unique()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        indata = adata.X.copy()
        genes = adata.var_names.copy()
        labels = _load_file_as_list(labels)
        if len(labels) != indata.shape[0]:
            raise ValueError(f"🛑 Length of training labels ({len(labels)}) does not match the number of input cells ({indata.shape[0]})")
    elif isinstance(X, str):
        raise ValueError("🛑 Invalid input. Supported types: .csv, .txt, .tsv, .tab, .mtx, .mtx.gz and .h5ad")
    else:
        logger.info("👀 The input training data is processed as an array-like object")
        indata = X.copy()
        if transpose:
            indata = indata.transpose()
        if isinstance(indata, pd.DataFrame):
            genes = indata.columns
            indata = indata.values
        else:
            if isinstance(indata, spmatrix):
                indata = indata.toarray()
            elif isinstance(indata, np.matrix):
                indata = np.array(indata)
            elif isinstance(indata, np.ndarray):
                indata = indata
            else:
                raise ValueError(f"🛑 Please provide a valid array-like object as input")
            if genes is None:
                raise Exception("🛑 Missing `genes`. Please provide this argument together with the input training data")
            genes = _load_file_as_list(genes)
            if len(genes) != indata.shape[1]:
                raise ValueError(f"🛑 The number of genes provided does not match the number of genes in the training data")
        if np.abs(np.expm1(indata[0]).sum()-10000) > 1:
            raise ValueError("🛑 Invalid expression matrix, expect log1p normalized expression to 10000 counts per cell")
        labels = _load_file_as_list(labels)
        if len(labels) != indata.shape[0]:
            raise ValueError(f"🛑 Length of training labels ({len(labels)}) does not match the number of input cells ({indata.shape[0]})")
    return indata, labels, genes

def _SGDClassifier(indata, labels,
                   alpha, max_iter, n_jobs,
                   mini_batch, batch_number, batch_size, epoch,
                   feature_selection, top_genes, **kwargs):
    """
    For internal use. Get the SGDClassifier.
    """
    classifier = SGDClassifier(loss = 'log', alpha = alpha, max_iter = max_iter, n_jobs = n_jobs, **kwargs)
    if not mini_batch:
        classifier.fit(indata, labels)
    else:
        return

def train(X = None,
          labels: Optional[Union[str, list, tuple, np.ndarray, pd.Series]] = None,
          genes: Optional[Union[str, list, tuple, np.ndarray, pd.Series]] = None,
          transpose_input: bool = False,
          #SGD param
          alpha: float = 0.0001, max_iter: int = 1000, n_jobs = None,
          #mini-batch
          mini_batch: bool = False, batch_number: int = 100, batch_size: int = 1000, epoch: int = 10,
          #feature selection
          feature_selection: bool = False, top_genes: int = 500,
          #other SGD param
          **kwargs
         ) -> Model:
    """
    coming soon...
    """
    logger.info("🍳 Preparing data before training")
    indata, labels, genes = _prepare_data(X, labels, genes, transpose_input)
    flag = indata.sum(axis = 0) == 0
    if flag.sum() > 0:
        logger.info(f"✂️ {flag.sum()} non-expressed genes are filtered out")
        indata = indata[:, ~flag]
        genes = genes[~flag]
    #scaler
    logger.info(f"⚖️ Scaling input data")
    scaler = StandardScaler()
    indata = scaler.fit_transform(indata)
    indata = np.clip(indata, a_min = None, a_max = 10)
    #classifier
    logger.info(f"🏋️ Training data")
    classifier = _SGDClassifier(indata = indata, labels = labels,
                                alpha = alpha, max_iter = max_iter, n_jobs = n_jobs,
                                mini_batch = mini_batch, batch_number = batch_number, batch_size = batch_size, epoch = epoch,
                                feature_selection = feature_selection, top_genes = top_genes, **kwargs)
