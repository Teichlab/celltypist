import os
import pickle
import requests
import numpy as np
import pandas as pd
from typing import List, Tuple
from celltypist import logger, defaults


data_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "data")

models_path = os.path.join(data_path, "models")


class Model():
    """Class that wraps SGDClassifier and its components."""
    def __init__(self, clf):
        self.classifier = clf
        if not hasattr(self.classifier, 'features'):
            with open(os.path.join(data_path, "genes.csv"), "rt") as f:
                self.classifier.features = np.array([gene.strip() for gene in f.readlines()])

    @staticmethod
    def load(model_file_path):
        if not os.path.exists(model_file_path):
            raise Exception(f"No such file: {model_file_path}")
        with open(model_file_path, "rb") as fh:
            try:
                return Model(pickle.load(fh)['Model'])
            except:
                raise Exception(f"Invalid model: {model_file_path}")

    @property
    def cells(self) -> List:
        """Get cells included in the model."""
        return list(self.classifier.classes_)
    
    # @property
    # def features(self) -> List:
    #     """Get genes included in the model."""
    #     with open(data_path, "rt") as f:
    #         return [gene.strip() for gene in f.readlines()]
    #     #return list(self.classifier.features)

    def predict_labels_and_prob(self, indata: np.ndarray) -> Tuple[List, List]:
        return self.classifier.predict(indata), self.classifier.predict_proba(indata)



def get_path_in_package(model: str) -> str:
    """Get a named model from the models included in the package."""
    return os.path.join(models_path, f"{model}.pkl")


def load(model: str = "") -> Model: #sklearn.linear_model.SGDClassifier:
    """Load model from package models, use model as named_model and fallback to using it as a path."""
    if not model:
        model = defaults.model
    if model in get_all_models():
        model = get_path_in_package(model)
    return Model.load(model)


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
    """Download models if there are non present in the models' directory"""
    if len([m for m in os.listdir(models_path) if m.endswith(".pkl")])==0:
        download_models()


def download_models() -> None:
    """Pull the model data from the server and download all the models"""
    url = 'https://celltypist.cog.sanger.ac.uk/models/models.json'
    logger.info("Retrieving model list from server")
    logger.info(f"Storing models in {models_path}")
    models_json = requests.get(url).json()
    for model in models_json["models"]:
        model_path = os.path.join(models_path, model["filename"])
        logger.info(f"downloading: {model['filename']}")
        try:
            with open(model_path, "wb") as f:
                f.write(requests.get(model["url"]).content)
        except Exception as exception:
            logger.error(f"{model['filename']} failed {exception}")


def update_models() -> None:
    """Update models by re-downloadig them."""
    download_models()

