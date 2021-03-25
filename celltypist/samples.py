import os

_samples_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "samples")


def _get_sample_data(filename: str) -> str:
    """Get the full path to the sample input data included in the package."""
    return os.path.join(_samples_path, filename)


def get_sample_csv() -> str:
    """
    Get the full path to the sample csv file included in the package.

    Returns
    ----------
    str
        A string of the full path to the sample csv file (`sample_cell_by_gene.csv`).
    """
    return _get_sample_data("sample_cell_by_gene.csv")
