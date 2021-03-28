import os
from typing import Optional, Literal, Union
import scanpy as sc
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from .models import Model
from . import logger
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
    decision_mat
        A :class:`~pandas.DataFrame` object returned from the celltyping process, showing the decision matrix.
    prob_mat
        A :class:`~pandas.DataFrame` object returned from the celltyping process, showing the probability matrix.
    adata
        An :class:`~scanpy.AnnData` object representing the input object.

    Attributes
    ----------
    predicted_labels
        Predicted labels including the individual prediction results and (if majority voting is done) majority voting results.
    decision_matrix
        Decision matrix with the decision score of each cell belonging to a given cell type.
    probability_matrix
        Probability matrix representing the probability each cell belongs to a given cell type (transformed from decision matrix by the sigmoid function).
    cell_count
        Number of input cells which undergo the prediction process.
    adata
        A Scanpy object representing the input data.
    """
    def __init__(self, labels: pd.DataFrame, decision_mat: pd.DataFrame, prob_mat: pd.DataFrame, adata: sc.AnnData):
        self.predicted_labels = labels
        self.decision_matrix = decision_mat
        self.probability_matrix = prob_mat
        self.adata = adata
        self.cell_count = labels.shape[0]

    def summary_frequency(self, by: Literal['predicted_labels', 'majority_voting'] = 'predicted_labels') -> pd.DataFrame:
        """
        Get the frequency of cells belonging to each cell type predicted by celltypist.

        Parameters
        ----------
        by
            Column name of :attr:`~celltypist.classifier.AnnotationResult.predicted_labels` specifying the prediction type which the summary is based on.
            Set to `majority_voting` if you want to summarize for the majority voting classifier.
            (Default: `predicted_labels`)

        Returns
        ----------
        :class:`~pandas.DataFrame`
            A :class:`~pandas.DataFrame` object with cell type frequencies.
        """
        unique, counts = np.unique(self.predicted_labels[by], return_counts=True)
        df = pd.DataFrame(list(zip(unique, counts)), columns=["celltype", "counts"])
        df.sort_values(['counts'], ascending=False, inplace=True)
        return df

    def to_adata(self, insert_labels: bool = True, insert_decision: bool = False, insert_probability: bool = False) -> sc.AnnData:
        """
        Insert the predicted labels, decision or probability matrix, and (if majority voting is done) majority voting results into the Scanpy object.

        Parameters
        ----------
        insert_labels
            Whether to insert the predicted cell type labels and (if majority voting is done) majority voting-based labels into the Scanpy object.
            (Default: `True`)
        insert_decision
            Whether to insert the decision matrix into the Scanpy object.
            (Default: `False`)
        insert_probability
            Whether to insert the probability matrix into the Scanpy object. This will override the decision matrix even when `insert_decision` is set to `True`.
            (Default: `False`)

        Returns
        ----------
        :class:`~scanpy.AnnData`
            Depending on whether majority voting is done, an :class:`~scanpy.AnnData` object with the following columns added to the observation metadata:
            1) **predicted_labels**, individual prediction outcome for each cell.
            2) **over_clustering**, over-clustering result for the cells.
            3) **majority_voting**, the cell type label assigned to each cell after the majority voting process.
            4) **name of each cell type**, which represents the decision scores (or probabilities if `insert_probability` is `True`) of a given cell type across cells.
        """
        if insert_labels:
            self.adata.obs[self.predicted_labels.columns] = self.predicted_labels
        if insert_probability:
            self.adata.obs[self.probability_matrix.columns] = self.probability_matrix
        elif insert_decision:
            self.adata.obs[self.decision_matrix.columns] = self.decision_matrix
        return self.adata

    def to_plots(self, folder: str, plot_probability: bool = False, format: str = 'pdf', prefix: str = '') -> None:
        """
        Plot the celltyping and (if majority voting is done) majority-voting results.

        Parameters
        ----------
        folder
            Path to a folder which stores the output figures.
        plot_probability
            Whether to also plot the decision score and probability distributions of each cell type across the test cells.
            If `True`, a number of figures will be generated (may take some time if the input data is large).
            (Default: `False`)
        format
            Format of output figures. Default to vector PDF files (note dots are still drawn with png backend).
            (Default: `pdf`)
        prefix
            Prefix for the output figures. Default to no prefix used.

        Returns
        ----------
        None
            Depending on whether majority voting is done and `plot_probability`, multiple UMAP plots showing the prediction and majority voting results in the `folder`:
            1) **predicted_labels**, individual prediction outcome for each cell overlaid onto the UMAP.
            2) **over_clustering**, over-clustering result of the cells overlaid onto the UMAP.
            3) **majority_voting**, the cell type label assigned to each cell after the majority voting process overlaid onto the UMAP.
            4) **name of each cell type**, which represents the decision scores and probabilities of a given cell type distributed across cells overlaid onto the UMAP.
        """
        if not os.path.isdir(folder):
            raise FileNotFoundError(f"🛑 Output folder {folder} does not exist. Please provide a valid folder")
        if 'X_umap' in self.adata.obsm:
            logger.info("👀 Detected existing UMAP coordinates, will plot the results accordingly")
        elif 'connectivities' in self.adata.obsp:
            logger.info("🧙 Generating UMAP coordinates based on the neighborhood graph")
            sc.tl.umap(self.adata)
        else:
            logger.info("🧙 Constructing the neighborhood graph and generating UMAP coordinates")
            adata = self.adata.copy()
            self.adata.obsm['X_pca'], self.adata.obsp['connectivities'], self.adata.obsp['distances'], self.adata.uns['neighbors'] = Classifier._construct_neighbor_graph(adata)
            sc.tl.umap(self.adata)
        logger.info("📈 Plotting the results")
        sc.settings.set_figure_params(figsize=[6.4, 6.4], format=format)
        self.adata.obs[self.predicted_labels.columns] = self.predicted_labels
        for column in self.predicted_labels:
            sc.pl.umap(self.adata, color = column, legend_loc = 'on data', show = False, legend_fontweight = 'normal', title = column.replace('_', ' '))
            plt.savefig(os.path.join(folder, prefix + column + '.' + format))
        if plot_probability:
            for column in self.probability_matrix:
                self.adata.obs['decision score'] = self.decision_matrix[column]
                self.adata.obs['probability'] = self.probability_matrix[column]
                sc.pl.umap(self.adata, color = ['decision score', 'probability'], show = False)
                plt.savefig(os.path.join(folder, prefix + column.replace('/','_') + '.' + format))
            self.adata.obs.drop(columns=['decision score', 'probability'], inplace=True)

    def to_table(self, folder: str, prefix: str = '', xlsx: bool = False) -> None:
        """
        Write out tables of predicted labels, decision matrix, and probability matrix.

        Parameters
        ----------
        folder
            Path to a folder which stores the output table/tables.
        prefix
            Prefix for the output table/tables. Default to no prefix used.
        xlsx
            Whether to merge output tables into a single Excel (.xlsx).
            (Default: `False`)

        Returns
        ----------
        None
            Depending on `xlsx`, return table(s) of predicted labels, decision matrix and probability matrix.
        """
        if not os.path.isdir(folder):
            raise FileNotFoundError(f"🛑 Output folder {folder} does not exist. Please provide a valid folder")
        if not xlsx:
            self.predicted_labels.to_csv(os.path.join(folder, f"{prefix}predicted_labels.csv"))
            self.decision_matrix.to_csv(os.path.join(folder, f"{prefix}decision_matrix.csv"))
            self.probability_matrix.to_csv(os.path.join(folder, f"{prefix}probability_matrix.csv"))
        else:
            with pd.ExcelWriter(os.path.join(folder, f"{prefix}annotation_result.xlsx")) as writer:
                self.predicted_labels.to_excel(writer, sheet_name="Predicted Labels")
                self.decision_matrix.to_excel(writer, sheet_name="Decision Matrix")
                self.probability_matrix.to_excel(writer, sheet_name="Probability Matrix")

    def __str__(self):
        return f"{self.cell_count} cells predicted into {len(np.unique(self.predicted_labels['predicted_labels']))} cell types"

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
        A Scanpy object which stores the log1p normalized expression data in `.X` or `.raw.X`.
    indata
        The expression matrix used for predictions stored in the log1p normalized format.
    indata_genes
        All the genes included in the input data.
    model
        A :class:`~celltypist.models.Model` object that wraps the SGDClassifier and the StandardScaler.
    """
    def __init__(self, filename: str, model: Model, transpose: bool = False, gene_file: Optional[str] = None, cell_file: Optional[str] = None):
        self.filename = filename
        logger.info(f"📁 Input file is '{self.filename}'")
        logger.info(f"⏳ Loading data")
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
            self.indata = self.adata.X.copy()
            self.indata_genes = self.adata.var_names.copy()
        elif self.filename.endswith('.h5ad'):
            self.adata = sc.read(self.filename)
            if self.adata.X.min() < 0:
                logger.info("👀 Detected scaled expression in the input data, will try the .raw attribute")
                try:
                    self.indata = self.adata.raw.X.copy()
                    self.indata_genes = self.adata.raw.var_names.copy()
                except Exception as e:
                    raise Exception(f"🛑 Fail to use the .raw attribute in the input object. {e}")
            else:
                self.indata = self.adata.X.copy()
                self.indata_genes = self.adata.var_names.copy()
            if np.abs(np.expm1(self.indata[0]).sum()-10000) > 1:
                raise ValueError("🛑 Invalid expression matrix, expect log1p normalized expression to 10000 counts per cell")
        else:
            raise ValueError("🛑 Invalid input file type. Supported types: .csv, .txt, .tsv, .tab, .mtx, .mtx.gz and .h5ad")

        logger.info(f"🔬 Input data has {self.indata.shape[0]} cells and {len(self.indata_genes)} genes")
        self.model = model

    def celltype(self) -> AnnotationResult:
        """
        Run celltyping jobs to predict cell types of input data.

        Returns
        ----------
        :class:`~celltypist.classifier.AnnotationResult`
            An :class:`~celltypist.classifier.AnnotationResult` object. Four important attributes within this class are:
            1) :attr:`~celltypist.classifier.AnnotationResult.predicted_labels`, predicted labels from celltypist.
            2) :attr:`~celltypist.classifier.AnnotationResult.decision_matrix`, decision matrix from celltypist.
            3) :attr:`~celltypist.classifier.AnnotationResult.probability_matrix`, probability matrix from celltypist.
            4) :attr:`~celltypist.classifier.AnnotationResult.adata`, Scanpy object representation of the input data.
        """
        logger.info(f"🔗 Matching reference genes in the model")
        k_x = np.isin(self.indata_genes, self.model.classifier.features)
        logger.info(f"🧬 {k_x.sum()} features used for prediction")
        k_x_idx = np.where(k_x)[0]
        self.indata = self.indata[:, k_x_idx]
        self.indata_genes = self.indata_genes[k_x_idx]
        lr_idx = pd.DataFrame(self.model.classifier.features, columns=['features']).reset_index().set_index('features').loc[self.indata_genes, 'index'].values

        logger.info(f"⚖️ Scaling input data")
        means_ = self.model.scaler.mean_[lr_idx]
        sds_ = self.model.scaler.scale_[lr_idx]
        self.indata = self.indata - means_
        self.indata = self.indata / sds_
        self.indata[self.indata > 10] = 10

        self.model.classifier.n_features_in_ = lr_idx.size
        self.model.classifier.features = self.model.classifier.features[lr_idx]
        self.model.classifier.coef_ = self.model.classifier.coef_[:, lr_idx]

        logger.info("🖋️ Predicting labels")
        decision_mat, prob_mat, lab = self.model.predict_labels_and_prob(self.indata)
        logger.info("✅ Prediction done!")

        cells = self.adata.obs_names
        return AnnotationResult(pd.DataFrame(lab, columns=['predicted_labels'], index=cells, dtype='category'), pd.DataFrame(decision_mat, columns=self.model.classifier.classes_, index=cells), pd.DataFrame(prob_mat, columns=self.model.classifier.classes_, index=cells), self.adata)

    @staticmethod
    def _construct_neighbor_graph(adata: sc.AnnData):
        """Construct a neighborhood graph. This function is for internal use."""
        if adata.X.min() < 0:
            adata = adata.raw.to_adata()
        sc.pp.filter_genes(adata, min_cells=5)
        sc.pp.highly_variable_genes(adata)
        adata = adata[:, adata.var.highly_variable]
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=50)
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
        return adata.obsm['X_pca'], adata.obsp['connectivities'], adata.obsp['distances'], adata.uns['neighbors']

    def over_cluster(self, resolution: Optional[float] = None) -> pd.Series:
        """
        Over-clustering input data with a canonical Scanpy pipeline. A neighborhood graph will be used (or constructed if not found) for the over-clustering.

        Parameters
        ----------
        resolution
            Resolution parameter for leiden clustering which controls the coarseness of the clustering.
            Default to 5, 10, 15 and 20 for datasets with cell numbers less than 5k, 20k, 40k and above, respectively.

        Returns
        ----------
        :class:`~pandas.Series`
            A :class:`~pandas.Series` object showing the over-clustering result.
        """
        if 'connectivities' not in self.adata.obsp:
            logger.info("👀 Can not detect a neighborhood graph, construct one before the over-clustering")
            adata = self.adata.copy()
            self.adata.obsm['X_pca'], self.adata.obsp['connectivities'], self.adata.obsp['distances'], self.adata.uns['neighbors'] = Classifier._construct_neighbor_graph(adata)
        else:
            logger.info("👀 Detected a neighborhood graph in the input object, will run over-clustering on the basis of it")
        if resolution is None:
            if self.adata.n_obs < 5000:
                resolution = 5
            elif self.adata.n_obs < 20000:
                resolution = 10
            elif self.adata.n_obs < 40000:
                resolution = 15
            else:
                resolution = 20
        logger.info(f"⛓️ Over-clustering input data with resolution set to {resolution}")
        sc.tl.leiden(self.adata, resolution=resolution, key_added='over_clustering')
        oc_column = self.adata.obs.over_clustering
        self.adata.obs.drop(columns=['over_clustering'], inplace=True)
        return oc_column

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
            An :class:`~celltypist.classifier.AnnotationResult` object. Four important attributes within this class are:
            1) :attr:`~celltypist.classifier.AnnotationResult.predicted_labels`, predicted labels from celltypist.
            2) :attr:`~celltypist.classifier.AnnotationResult.decision_matrix`, decision matrix from celltypist.
            3) :attr:`~celltypist.classifier.AnnotationResult.probability_matrix`, probability matrix from celltypist.
            4) :attr:`~celltypist.classifier.AnnotationResult.adata`, Scanpy object representation of the input data.
        """
        if isinstance(over_clustering, list):
            over_clustering = np.array(over_clustering)
        logger.info("🗳️ Majority voting the predictions")
        votes = pd.crosstab(predictions.predicted_labels['predicted_labels'], over_clustering)
        majority = votes.idxmax()[over_clustering].reset_index()
        majority.index = predictions.predicted_labels.index
        majority.columns = ['over_clustering', 'majority_voting']
        majority['majority_voting'] = majority['majority_voting'].astype('category')
        predictions.predicted_labels = predictions.predicted_labels.join(majority)
        logger.info("✅ Majority voting done!")
        return predictions
