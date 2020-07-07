import os
import math
from typing import Tuple
import numpy as np
import pandas as pd
# parallelisation
from joblib import Parallel, delayed
# hide warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)


class AnnotationResult():
    """Class that represents the result of a celltyping annotation process."""
    def __init__(self, labels: np.ndarray, prob_matrix: np.ndarray, model):
        self.predicted_labels = pd.DataFrame(labels)
        self.probability_matrix = pd.DataFrame(prob_matrix, columns=model.classes_)
        self.model = model
        self.cell_count = labels.shape
        self.summary = self.summarize_labels(labels)

    def summarize_labels(self, labels: np.ndarray) -> pd.DataFrame:
        """Get a summary of the cells per label obtained in the annotation process."""
        unique, counts = np.unique(labels, return_counts=True)
        df = pd.DataFrame(list(zip(unique, counts)), columns=["celltype", "counts"])
        df.sort_values(['counts'], ascending=False, inplace=True)
        return df

    def write_excel(self, filename: str):
        """Write excel file with both the predicted labels and the probability matrix."""
        filename, _ = os.path.splitext(filename)
        with pd.ExcelWriter(f"{filename}.xlsx") as writer:
            self.predicted_labels.to_excel(writer, sheet_name="Predicted Labels")
            self.probability_matrix.to_excel(writer, sheet_name="Probability Matrix")


class Classifier():
    """Class that wraps the cell typing process."""
    def __init__(self, filename: str, model, chunk_size: int, cpus: int, quiet: bool):
        self.filename = filename
        self.chunk_size = chunk_size
        self.cpus = cpus
        self.model = model
        with open(self.filename) as fh:
            self.cell_count = sum(1 for line in fh)
        self.chunk_iterator = range(math.ceil(self.cell_count/self.chunk_size))
        self.quiet = quiet

    def process_chunk(self, start_at: int) -> Tuple[np.ndarray, np.ndarray]:
        """Process a chunk of the input file starting at the offset position."""
        X_test = np.log1p(pd.read_csv(self.filename, skiprows=start_at, nrows=self.chunk_size, header=None, index_col=0).values)
        return self.model.predict(X_test), self.model.predict_proba(X_test)

    def celltype(self) -> AnnotationResult:
        """Run celltyping jobs to get results."""
        result = Parallel(n_jobs=self.cpus, verbose=10 if not self.quiet else 0)(
            delayed(self.process_chunk)(start_at=i*self.chunk_size+1) for i in self.chunk_iterator)
        lab_mat = np.hstack([result[i][0] for i in range(len(result))])
        prob_mat = np.vstack([result[i][1] for i in range(len(result))])
        return AnnotationResult(lab_mat, prob_mat, self.model)

    def print_config(self):
        """Show current configuration values for this clasifier."""
        (f"filename={self.filename}. cpus={self.cpus}. chunk_size={self.chunk_size}")
