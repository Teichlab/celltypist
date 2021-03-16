import os
from typing import Optional, Literal, Union
import scanpy as sc
import numpy as np
import pandas as pd
from .models import Model
from . import logger
# parallelisation
#from joblib import Parallel, delayed
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)


class AnnotationResult():
    """
    Class that represents the result of a celltyping annotation process.

    Parameters
    ----------
    labels
        A :class:`~pandas.DataFrame` object returned from the celltyping process, showing the predicted labels.
    prob
        A :class:`~pandas.DataFrame` object returned from the celltyping process, showing the probability matrix.

    Attributes
    ----------
    predicted_labels
        Predicted labels including the individual prediction results and (if majority voting is done) majority voting results.
    probability_table
        Probability matrix representing the probability each cell belongs to a given cell type.
    cell_count
        Number of input cells which are predicted by celltypist.
    """
    def __init__(self, labels: pd.DataFrame, prob: pd.DataFrame):
        self.predicted_labels = labels
        self.probability_table = prob
        self.cell_count = labels.shape[0]

    def summary_frequency(self, by: Literal['predicted labels', 'majority voting'] = 'predicted labels') -> pd.DataFrame:
        """
        Get the frequency of cells belonging to each cell type predicted by celltypist.

        Parameters
        ----------
        by
            Column name of :attr:`~celltypist.classifier.AnnotationResult.predicted_labels` specifying the prediction type which the summary is based on.
            Set to `majority voting` if you want to summarize for the majority voting classifier.
            (Default: `predicted labels`)

        Returns
        ----------
        :class:`~pandas.DataFrame`
            A :class:`~pandas.DataFrame` object with cell type frequencies.
        """
        unique, counts = np.unique(self.predicted_labels[by], return_counts=True)
        df = pd.DataFrame(list(zip(unique, counts)), columns=["celltype", "counts"])
        df.sort_values(['counts'], ascending=False, inplace=True)
        return df

    def write_excel(self, filename: str) -> None:
        """
        Write excel file with both the predicted labels and the probability matrix.

        Parameters
        ----------
        filename
            Excel file (.xlsx) to store the predicted cell types and probability matrix.

        Returns
        ----------
        None
            xlsx table containing two sheets of predicted labels and probability matrix, respectively.
        """
        filename, _ = os.path.splitext(filename)
        with pd.ExcelWriter(f"{filename}.xlsx") as writer:
            self.predicted_labels.to_excel(writer, sheet_name="Predicted Labels")
            self.probability_table.to_excel(writer, sheet_name="Probability Matrix")

    def __str__(self):
        return f"{self.cell_count} cells predicted into {len(np.unique(self.predicted_labels['predicted labels']))} cell types"

class Classifier():
    """
    Class that wraps the celltyping and majority voting processes.

    Parameters
    ----------
    filename
        Path to the input count matrix (supported types are csv, txt, tsv, tab and mtx) or Scanpy object (h5ad).
        If it's the former, a cell-by-gene format is desirable (see `transpose` for more information).
        Genes should be gene symbols. Non-expressed genes are preferred to be provided as well.
    model
        A :class:`~celltypist.models.Model` object that wraps the SGDClassifier and the StandardScaler.
    transpose
        Whether to transpose the input matrix. Set to `True` if `filename` is provided in a gene-by-cell format.
        (Default: `False`)
    gene_file
        Path to the file which stores each gene per line corresponding to the genes used in the provided mtx file.
        Ignored if `filename` is not provided in the mtx format.
    cell_file
        Path to the file which stores each cell per line corresponding to the cells used in the provided mtx file.
        Ignored if `filename` is not provided in the mtx format.

    Attributes
    ----------
    filename
        Path to the input dataset.
    adata
        A Scanpy object which stores the log1p normalized expression data.
    indata
        The expression matrix used for predictions stored in the log1p normalized format.
    indata_genes
        All the genes included in the input data.
    model
        A :class:`~celltypist.models.Model` object that wraps the SGDClassifier and the StandardScaler.
    """
    def __init__(self, filename: str, model: Model, transpose: bool = False, gene_file: Optional[str] = None, cell_file: Optional[str] = None): #, chunk_size: int, cpus: int, quiet: bool):
        self.filename = filename
        logger.info(f"📁 Input file is '{self.filename}'")
        logger.info(f"⏳ Loading data...")
        if self.filename.endswith(('.csv', '.txt', '.tsv', '.tab', '.mtx', '.mtx.gz')):
            self.adata = sc.read(self.filename)
            if transpose:
                self.adata = self.adata.transpose()
            if self.filename.endswith(('.mtx', '.mtx.gz')):
                if (gene_file is None) or (cell_file is None):
                    raise FileNotFoundError("🛑 Missing `gene_file` and/or `cell_file`. Please provide both arguments together with the input mtx file")
                genes_mtx = pd.read_csv(gene_file, header=None)[0].values
                cells_mtx = pd.read_csv(cell_file, header=None)[0].values
                if len(genes_mtx) != self.adata.n_vars:
                    raise ValueError(f"🛑 The number of genes in {gene_file} does not match the number of genes in {self.filename}")
                if len(cells_mtx) != self.adata.n_obs:
                    raise ValueError(f"🛑 The number of cells in {cell_file} does not match the number of cells in {self.filename}")
                self.adata.var_names = genes_mtx
                self.adata.obs_names = cells_mtx
            self.adata.var_names_make_unique()
            sc.pp.normalize_total(self.adata, target_sum=1e4)
            sc.pp.log1p(self.adata)
        elif self.filename.endswith('.h5ad'):
            self.adata = sc.read(self.filename)
            if self.adata.X.min() < 0:
                logger.info("👀 Detect scaled expression in the input data, will try the .raw attribute...")
                try:
                    self.adata = self.adata.raw.to_adata()
                except Exception:
                    raise Exception("🛑 Fail to use the .raw attribute in the input object")
            if np.abs(np.expm1(self.adata.X[0]).sum()-10000) > 1:
                raise ValueError("🛑 Invalid expression matrix, expect log1p normalized expression to 10000 counts per cell")
        else:
            raise ValueError("🛑 Invalid input file type. Supported types: .csv, .txt, .tsv, .tab, .mtx, .mtx.gz and .h5ad")
        self.indata = self.adata.X.copy()
        self.indata_genes = self.adata.var_names.copy()

        logger.info(f"🔬 Input data has {self.indata.shape[0]} cells and {len(self.indata_genes)} genes")
        # self.chunk_size = chunk_size
        # self.cpus = cpus
        self.model = model
        #with open(self.filename) as fh:
        #    self.cell_count = sum(1 for line in fh)
        #self.chunk_iterator = range(math.ceil(self.cell_count/self.chunk_size))
        #self.quiet = quiet

    #def process_chunk(self, start_at: int) -> None: #-> Tuple[np.ndarray, np.ndarray]:
        #"""Process a chunk of the input file starting at the offset position."""
        #X_test = np.log1p(pd.read_csv(self.filename, skiprows=start_at, nrows=self.chunk_size, header=None, index_col=0).values)
        #return self.model.predict_labels_and_prob(X_test)
        #pass

    def celltype(self) -> AnnotationResult:
        """
        Run celltyping jobs to predict cell types of input data.

        Returns
        ----------
        :class:`~celltypist.classifier.AnnotationResult`
            An :class:`~celltypist.classifier.AnnotationResult` object. Two important attributes within this class are:
            1) :attr:`~celltypist.classifier.AnnotationResult.predicted_labels`, predicted labels from celltypist.
            2) :attr:`~celltypist.classifier.AnnotationResult.probability_table`, probability matrix from celltypist.
        """
        #result = Parallel(n_jobs=self.cpus, verbose=10 if not self.quiet else 0)(
        #    delayed(self.process_chunk)(start_at=i*self.chunk_size+1) for i in self.chunk_iterator)
        #lab_mat = np.hstack([result[i][0] for i in range(len(result))])
        #prob_mat = np.vstack([result[i][1] for i in range(len(result))])

        logger.info(f"🧙 Matching reference genes")
        k_x = np.isin(self.indata_genes, self.model.classifier.features)
        logger.info(f"🧩 {k_x.sum()} features used for prediction")
        k_x_idx = np.where(k_x)[0]
        self.indata = self.indata[:, k_x_idx]
        self.indata_genes = self.indata_genes[k_x_idx]
        lr_idx = pd.DataFrame(self.model.classifier.features, columns=['features']).reset_index().set_index('features').loc[self.indata_genes, 'index'].values

        logger.info(f"🧙 Scaling input data")
        means_ = self.model.scaler.mean_[lr_idx]
        sds_ = self.model.scaler.scale_[lr_idx]
        self.indata = self.indata - means_
        self.indata = self.indata / sds_

        self.model.classifier.n_features_in_ = lr_idx.size
        self.model.classifier.features = self.model.classifier.features[lr_idx]
        self.model.classifier.coef_ = self.model.classifier.coef_[:, lr_idx]

        logger.info("🖋️ Predicting labels")
        lab_mat, prob_mat = self.model.predict_labels_and_prob(self.indata)
        # print(results)
        # # lab_mat = np.hstack(results)
        # # prob_mat = np.vstack(results)
        logger.info("✅ Prediction done!")

        cells = self.adata.obs_names
        return AnnotationResult(pd.DataFrame(lab_mat, columns=['predicted labels'], index=cells), pd.DataFrame(prob_mat, columns=self.model.classifier.classes_, index=cells))

    def over_cluster(self, resolution: Optional[int] = None) -> pd.Series:
        """
        Over-clustering input data with a canonical scanpy pipeline.

        Parameters
        ----------
        resolution
            resolution parameter for leiden clustering which controls the coarseness of the clustering.
            Default to 5, 10, 15 and 20 for datasets with cell numbers less than 5k, 20k, 40k and above, respectively.

        Returns
        ----------
        :class:`~pandas.Series`
            A :class:`~pandas.Series` object showing the over-clustering result.
        """
        if resolution is None:
            if self.adata.shape[0] < 5000:
                resolution = 5
            elif self.adata.shape[0] < 20000:
                resolution = 10
            elif self.adata.shape[0] < 40000:
                resolution = 15
            else:
                resolution = 20
        logger.info(f"🧙 Over-clustering input data with resolution set to {resolution}")
        sc.pp.filter_genes(self.adata, min_cells=1)
        sc.pp.highly_variable_genes(self.adata)
        sc.pp.scale(self.adata, max_value=10)
        sc.tl.pca(self.adata, n_comps=50)
        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=50)
        sc.tl.leiden(self.adata, resolution=resolution, key_added='over_clustering')
        return self.adata.obs['over_clustering']

    @staticmethod
    def majority_vote(predictions: AnnotationResult, over_clustering: Union[list, np.ndarray, pd.Series]) -> AnnotationResult:
        """
        Majority vote the celltypist predictions using the result from the over-clustering.

        Parameters
        ----------
        predictions
            An :class:`~celltypist.classifier.AnnotationResult` object containing the :attr:`~celltypist.classifier.AnnotationResult.predicted_labels`.
        over_clustering
            A list, numpy array or pandas series containing the over-clustering information.

        Returns
        ----------
        :class:`~celltypist.classifier.AnnotationResult`
            An :class:`~celltypist.classifier.AnnotationResult` object. Two important attributes within this class are:
            1) :attr:`~celltypist.classifier.AnnotationResult.predicted_labels`, predicted labels from celltypist.
            2) :attr:`~celltypist.classifier.AnnotationResult.probability_table`, probability matrix from celltypist.
        """
        if isinstance(over_clustering, list):
            over_clustering = np.array(over_clustering)
        logger.info("🧙 Majority voting")
        votes = pd.crosstab(predictions.predicted_labels['predicted labels'], over_clustering)
        majority = votes.idxmax()[over_clustering].reset_index()
        majority.index = predictions.predicted_labels.index
        majority.columns = ['over clustering', 'majority voting']
        predictions.predicted_labels = predictions.predicted_labels.join(majority)
        logger.info("✅ Majority voting done!")
        return predictions

    # def print_config(self):
    #     """Show current configuration values for this clasifier."""
    #     (f"filename={self.filename}. cpus={self.cpus}. chunk_size={self.chunk_size}")
