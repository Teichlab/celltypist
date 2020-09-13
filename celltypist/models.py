import os
import pickle
from typing import List
import requests
from celltypist import logger

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
        return pickle.load(fh)['Model']


def get_all_models() -> List[str]:
    """Get a List of all the avaiable models included in the package."""
    download_if_required()
    avaiable_models = []
    for model_filename in os.listdir(models_path):
        if model_filename.endswith(".pkl"):
            model_name = os.path.basename(model_filename)[:-4]
            avaiable_models.append(model_name)        
    return avaiable_models


def download_if_required() -> None:
      if len([m for m in os.listdir(models_path) if m.endswith(".pkl")])==0:
          download_models()

def download_models() -> None:
    url = 'https://celltypist.cog.sanger.ac.uk/models/models.json'
    logger.info("Retrieving model list from server")
    models_json = requests.get(url).json()
    for model in models_json["models"]:
        logger.info(f"Downloading: {model['filename']}")
        with open(os.path.join(models_path, model["filename"]), "wb") as f:
            f.write(requests.get(model["url"]).content)

def update_models() -> None:
    download_models()