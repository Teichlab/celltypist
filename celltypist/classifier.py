import os
import math
import anndata
from typing import Tuple
import numpy as np
import pandas as pd
from celltypist.models import Model
from celltypist import models, helpers, logger
# parallelisation
from joblib import Parallel, delayed
# hide warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
import scanpy as sc


class AnnotationResult():
    """Class that represents the result of a celltyping annotation process."""
    def __init__(self, labels: np.ndarray, prob_matrix: np.ndarray, model, indata):
        self.predicted_labels = labels
        self.probability_matrix = prob_matrix
        self.cell_count = labels.shape
        self.model_celltypes = model.classifier.classes_
        self.indata = indata

    def predicted_labels_as_df(self) -> pd.DataFrame:
        return pd.DataFrame(self.predicted_labels)

    def probability_matrix_as_df(self) -> pd.DataFrame:
        return pd.DataFrame(self.probability_matrix, columns=self.model_celltypes)

    def summary_as_df(self) -> pd.DataFrame:
        """Get a summary of the cells per label obtained in the annotation process."""
        unique, counts = np.unique(self.predicted_labels, return_counts=True)
        df = pd.DataFrame(list(zip(unique, counts)), columns=["celltype", "counts"])
        df.sort_values(['counts'], ascending=False, inplace=True)
        return df

    def write_excel(self, filename: str):
        """Write excel file with both the predicted labels and the probability matrix."""
        filename, _ = os.path.splitext(filename)
        with pd.ExcelWriter(f"{filename}.xlsx") as writer:
            self.predicted_labels_as_df().to_excel(writer, sheet_name="Predicted Labels")
            self.probability_matrix_as_df().to_excel(writer, sheet_name="Probability Matrix")
 
    def __str__(self):
        return f"cell count: {self.cell_count}\npredicted labels: {self.predicted_labels}\nmodel celltypist{self.model_celltypes}"

class Classifier():
    """Class that wraps the cell typing process."""
    def __init__(self, filename: str, model: Model, transpose = False): #, chunk_size: int, cpus: int, quiet: bool):
        self.filename = filename
        logger.info(f"📁 Input file is '{self.filename}'")
        logger.info(f"⏳ Loading data...")
        if self.filename.endswith('.csv'):
            self.file_type = "csv"
            self.adata = sc.read(self.filename)
            if transpose:
                self.adata = self.adata.transpose()
            sc.pp.normalize_total(self.adata)
            sc.pp.log1p(self.adata)
        elif helpers.is_h5ad(self.filename):
            self.file_type = "h5ad"
            self.adata = sc.read(self.filename)
            if self.adata.X.min() < 0:
                raise Exception("🛑 Detect scaled expression while expect log1p normalized expression (sc.pp.normalize_total followed by sc.pp.log1p)")
        else:
            raise Exception("🛑 Invlaid input file type. Supported types: .csv and .h5ad")
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

    def process_chunk(self, start_at: int) -> None: #-> Tuple[np.ndarray, np.ndarray]:
        """Process a chunk of the input file starting at the offset position."""
        #X_test = np.log1p(pd.read_csv(self.filename, skiprows=start_at, nrows=self.chunk_size, header=None, index_col=0).values)
        #return self.model.predict_labels_and_prob(X_test)
        pass

    def celltype(self) -> AnnotationResult:
        """Run celltyping jobs to get results."""
        #result = Parallel(n_jobs=self.cpus, verbose=10 if not self.quiet else 0)(
        #    delayed(self.process_chunk)(start_at=i*self.chunk_size+1) for i in self.chunk_iterator)
        #lab_mat = np.hstack([result[i][0] for i in range(len(result))])
        #prob_mat = np.vstack([result[i][1] for i in range(len(result))])
        
        logger.info(f"🧙 Gene reference matching")
        ############################################################
        #def gene_reference_matching(self, input_data, input_genes):
        ############################################################    
        # features is the gene names of the dataset to be annotated
        k_x = np.isin(self.indata_genes, list(self.model.classifier.features))
        logger.info(f"🧩 {k_x.sum()} features used for prediction")
        k_x_idx = np.where(k_x)[0]
        self.indata = self.indata[:, k_x_idx]
        self.indata_genes = self.indata_genes[k_x]

        lr_idx = pd.DataFrame(self.model.classifier.features, columns=['features']).reset_index().set_index('features').loc[self.indata_genes, 'index'].values

        means_ = self.model.scaler.mean_[lr_idx]
        sds_ = self.model.scaler.scale_[lr_idx]
        self.indata = self.indata - means_
        self.indata = self.indata / sds_
        #self.indata[:, self.indata.std(axis=0).A1 == 0] = self.indata.min()

        self.model.classifier.n_features_in_ = lr_idx.size
        self.model.classifier.features = self.model.classifier.features[lr_idx]
        self.model.classifier.coef_ = self.model.classifier.coef_[:, lr_idx]

        logger.info("🖋️ Predicting labels")
        lab_mat, prob_mat = self.model.predict_labels_and_prob(self.indata)
        # print(results)
        # # lab_mat = np.hstack(results)
        # # prob_mat = np.vstack(results)
        logger.info("✅ Done!")

        return AnnotationResult(lab_mat, prob_mat, self.model, self.indata)

    # def print_config(self):
    #     """Show current configuration values for this clasifier."""
    #     (f"filename={self.filename}. cpus={self.cpus}. chunk_size={self.chunk_size}")
