import os

samples_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "samples")


def get_sample_data(filename: str) -> str:
    """Get path to sample input data included in the pacakge."""
    return os.path.join(samples_path, filename)


def get_sample_csv() -> str:
    """Get path for the sample CSV file included in the package."""
    return get_sample_data("sample_cell_by_gene.csv")
