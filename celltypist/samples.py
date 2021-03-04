import os

samples_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "samples")


def get_sample_data(filename: str) -> str:
    """Get the full path to the sample input data included in the package."""
    return os.path.join(samples_path, filename)


def get_sample_csv() -> str:
    """
    Get the full path to the sample csv file included in the package.

    Returns
    ----------
    A string of the full path to the sample csv file (`sample_cell_by_gene.csv`).
    """
    return get_sample_data("sample_cell_by_gene.csv")
