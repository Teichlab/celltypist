import os
import json
import pickle
import requests
import numpy as np
import pandas as pd
from typing import List, Tuple
from . import logger


data_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "data")

models_path = os.path.join(data_path, "models")


class Model():
    """Class that wraps SGDClassifier and its components."""
    def __init__(self, clf, scaler):
        self.classifier = clf
        if not hasattr(self.classifier, 'features'):
            with open(os.path.join(data_path, "genes.csv"), "rt") as f:
                self.classifier.features = np.array([gene.strip() for gene in f.readlines()])
        self.scaler = scaler

    @staticmethod
    def load(model_file_path):
        if not os.path.exists(model_file_path):
            raise FileNotFoundError(f"🛑 No such file: {model_file_path}")
        with open(model_file_path, "rb") as fh:
            try:
                pkl_obj = pickle.load(fh)
                return Model(pkl_obj['Model'], pkl_obj['Scaler_'])
            except Exception as exception:
                raise Exception(f"🛑 Invalid model: {model_file_path}. {exception}")

    @property
    def cell_types(self) -> List:
        """Get cell types included in the model."""
        return self.classifier.classes_

    @property
    def features(self) -> List:
        """Get genes included in the model."""
        return self.classifier.features

    def predict_labels_and_prob(self, indata: np.ndarray) -> Tuple[List, List]:
        return self.classifier.predict(indata), self.classifier.predict_proba(indata)


def get_path_in_package(model: str) -> str:
    """
    Get the full path to the model specified.

    Parameters
    ----------
    model
        Model name as a string.
        To see all available models and their descriptions, use `celltypist.models.models_description()`

    Returns
    ----------
    A string of the full path to the desired model.
    """
    return os.path.join(models_path, f"{model}")


def load(model: str = "") -> Model:
    """
    Load the desired model.

    Parameters
    ----------
    model
        Model name specifying the model you want to load. Default to `Immune_All_Low.pkl` if not provided.
        To see all available models and their descriptions, use `celltypist.models.models_description()`.

    Returns
    ----------
    A `~celltypist.models.Model` object.
    """
    if not model:
        model = get_default_model()
    if model in get_all_models():
        model = get_path_in_package(model)
    return Model.load(model)


def get_default_model() -> str:
    """
    Get the default model name.

    Returns
    ----------
    A string showing the default model name (should be `Immune_All_Low.pkl`).
    """
    models_json = get_models_index()
    default_model = [m["filename"] for m in models_json["models"] if ("default" in m and m["default"])]
    if not default_model:
        first_model = models_json["models"][0]["filename"]
        logger.warn(f"👀 No model marked as 'default', using {first_model}")
        return first_model
    if len(default_model) > 1:
        logger.warn(f"👀 More than one model marked as 'default', using {default_model[0]}")
    return default_model[0]


def get_all_models() -> List[str]:
    """
    Get a list of all the available models included in the package.

    Returns
    ----------
    A list of available models.
    """
    download_if_required()
    avaiable_models = []
    for model_filename in os.listdir(models_path):
        if model_filename.endswith(".pkl"):
            model_name = os.path.basename(model_filename)
            avaiable_models.append(model_name)
    return avaiable_models


def download_if_required() -> None:
    """Download models if there are none present in the package models directory."""
    if len([m for m in os.listdir(models_path) if m.endswith(".pkl")]) == 0:
        logger.info(f"🔎 No available models. Downloading...")
        download_models()


def get_models_index(force_update: bool=False) -> dict:
    """
    Get the model json object containing the model list.

    Parameters
    ----------
    force_update
        If set to `True`, will download the latest model json file from the remote.
        (Default: `False`)

    Returns
    ----------
    A dict object converted from the json file.
    """
    models_json_path = os.path.join(models_path, "models.json")
    if not os.path.exists(models_json_path) or force_update:
        download_model_index()
    with open(models_json_path) as f:
        return json.load(f)


def download_model_index(only_model: bool = True) -> None:
    """
    Download the models.json file from the remote server.

    Parameters
    ----------
    only_model
        If set to `False`, will also download the models in addition to the json file.
        (Default: `True`)
    """
    url = 'https://celltypist.cog.sanger.ac.uk/models/models.json'
    logger.info(f"📜 Retrieving model list from server {url}")
    with open(os.path.join(models_path, "models.json"), "wb") as f:
        f.write(requests.get(url).content)
    model_count = len(requests.get(url).json()["models"])
    logger.info(f"📚 Total models in list: {model_count}")
    if not only_model:
        download_models()

def download_models(force_update: bool=False) -> None:
    """
    Download all the available models.

    Parameters
    ----------
    force_update
        Whether to fetch a latest JSON index for downloading all available models.
        Set to `True` if you want to parallel the latest celltypist model releases.
        (Default: `False`)
    """
    models_json = get_models_index(force_update)
    model_count = len(models_json["models"])
    logger.info(f"📂 Storing models in {models_path}")
    for idx,model in enumerate(models_json["models"]):
        model_path = os.path.join(models_path, model["filename"])
        logger.info(f"💾 Downloading model [{idx+1}/{model_count}]: {model['filename']}")
        try:
            with open(model_path, "wb") as f:
                f.write(requests.get(model["url"]).content)
        except Exception as exception:
            logger.error(f"🛑 {model['filename']} failed {exception}")


def update_models() -> None:
    """Update models by re-downloading the latest model versions."""
    download_models(force_update=True)

def models_description() -> pd.DataFrame:
    """
    Get the descriptions of all available models.

    Returns
    ----------
    A `~pandas.DataFrame` object with model descriptions.
    """
    models_json = get_models_index()
    models = models_json["models"]
    filenames = [model['filename'] for model in models]
    descriptions = [model['details'] for model in models]
    return pd.DataFrame({'model': filenames, 'description': descriptions})
