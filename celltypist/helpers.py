import csv
import numpy as np
import pandas as pd

def prepare():
    """Prepare dataset to be annotated."""
    pass


def validate():
    """Validte dataset can be processed."""
    pass


# def get_gene_names(input_filename):
#     """Get gene names from the file."""
#     df = pd.read_csv(input_filename, header=None, index_col=0, nrows=1)
#     return df.to_numpy()


def is_h5ad(input_file: str) -> bool:
    with open(input_file, "rb") as f:
        return str(f.read(5)[1:4]) == "b'HDF'"

#def is_csv(input_file: str) -> bool:
#    try:
#        df = pd.read_csv(input_file, header=None, index_col=0, nrows=1)
#        return True
#    except:
#        return False
