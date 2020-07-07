import os
import pickle
from typing import List

models_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "data", "models")


def get_path_in_package(model: str) -> str:
    """Get a named model from the models included in the package."""
    return os.path.join(models_path, f"{model}.pkl")


def load(model: str = "default") -> str:
    """Load model from package models, use model as named_model and fallback to using it as a path."""
    if model in get_all_models():
        model = get_path_in_package(model)
    with open(model, "rb") as fh:
        return pickle.load(fh)


def get_all_models() -> List[str]:
    """Get a List of all the avaiable models included in the package."""
    avaiable_models = []
    for model_filename in os.listdir(models_path):
        if model_filename.endswith(".pkl"):
            model_name = os.path.basename(model_filename)[:-4]
            avaiable_models.append(model_name)
    return avaiable_models
