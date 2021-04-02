import os
import pathlib
import json
import pickle
import requests
import numpy as np
import pandas as pd
from typing import Optional
from scipy.special import expit
from . import logger

#create ~/.celltypist and subdirs
celltypist_path = os.path.join(str(pathlib.Path.home()), '.celltypist')
pathlib.Path(celltypist_path).mkdir(parents=True, exist_ok=True)
data_path = os.path.join(celltypist_path, "data")
models_path = os.path.join(data_path, "models")
pathlib.Path(models_path).mkdir(parents=True, exist_ok=True)


class Model():
    """
    Class that wraps the SGDClassifier and the StandardScaler.

    Parameters
    ----------
    clf
        A SGDClassifier incorporated in the loaded model.
    scaler
        A StandardScaler incorporated in the loaded model.

    Attributes
    ----------
    classifier
        The SGDClassifier incorporated in the loaded model.
    scaler
        The StandardScaler incorporated in the loaded model.
    """
    def __init__(self, clf, scaler):
        self.classifier = clf
        self.scaler = scaler

    @staticmethod
    def load(model_file_path: str):
        """Load the desired model."""
        if not os.path.isfile(model_file_path):
            raise FileNotFoundError(f"🛑 No such file: {model_file_path}")
        with open(model_file_path, "rb") as fh:
            try:
                pkl_obj = pickle.load(fh)
                return Model(pkl_obj['Model'], pkl_obj['Scaler_'])
            except Exception as exception:
                raise Exception(f"🛑 Invalid model: {model_file_path}. {exception}")

    @property
    def cell_types(self) -> np.ndarray:
        """Get cell types included in the model."""
        return self.classifier.classes_

    @property
    def features(self) -> np.ndarray:
        """Get genes included in the model."""
        return self.classifier.features

    def predict_labels_and_prob(self, indata) -> tuple:
        """Get the decision matrix, probability matrix, and predicted cell types for the input data."""
        scores = self.classifier.decision_function(indata)
        probs = expit(scores)
        return scores, probs, self.classifier.classes_[scores.argmax(axis=1)]


def get_model_path(file: str) -> str:
    """
    Get the full path to a file in the `models` folder.

    Parameters
    ----------
    file
        File name as a string.
        To see all available models and their descriptions, use :func:`~celltypist.models.models_description()`.

    Returns
    ----------
    str
        A string of the full path to the desired file.
    """
    return os.path.join(models_path, f"{file}")


def load(model: Optional[str] = None) -> Model:
    """
    Load the desired model.

    Parameters
    ----------
    model
        Model name specifying the model you want to load. Default to `Immune_All_Low.pkl` if not provided.
        To see all available models and their descriptions, use :func:`~celltypist.models.models_description()`.

    Returns
    ----------
    :class:`~celltypist.models.Model`
        A :class:`~celltypist.models.Model` object.
    """
    if model is None:
        model = get_default_model()
    if model in get_all_models():
        model = get_model_path(model)
    return Model.load(model)


def get_default_model() -> str:
    """
    Get the default model name.

    Returns
    ----------
    str
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


def get_all_models() -> list:
    """
    Get a list of all the available models.

    Returns
    ----------
    list
        A list of available models.
    """
    download_if_required()
    available_models = []
    for model_filename in os.listdir(models_path):
        if model_filename.endswith(".pkl"):
            model_name = os.path.basename(model_filename)
            available_models.append(model_name)
    return available_models


def download_if_required() -> None:
    """Download models if there are none present in the `models` directory."""
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
    dict
        A dict object converted from the model json file.
    """
    models_json_path = get_model_path("models.json")
    if not os.path.exists(models_json_path) or force_update:
        download_model_index()
    with open(models_json_path) as f:
        return json.load(f)


def download_model_index(only_model: bool = True) -> None:
    """
    Download the `models.json` file from the remote server.

    Parameters
    ----------
    only_model
        If set to `False`, will also download the models in addition to the json file.
        (Default: `True`)
    """
    url = 'https://celltypist.cog.sanger.ac.uk/models/models.json'
    logger.info(f"📜 Retrieving model list from server {url}")
    with open(get_model_path("models.json"), "wb") as f:
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
        model_path = get_model_path(model["filename"])
        if os.path.exists(model_path) and not force_update:
            logger.info(f"⏩ Skipping [{idx+1}/{model_count}]: {model['filename']} (file exists)")
            continue
        logger.info(f"💾 Downloading model [{idx+1}/{model_count}]: {model['filename']}")
        try:
            with open(model_path, "wb") as f:
                f.write(requests.get(model["url"]).content)
        except Exception as exception:
            logger.error(f"🛑 {model['filename']} failed {exception}")


def models_description() -> pd.DataFrame:
    """
    Get the descriptions of all available models.

    Returns
    ----------
    :class:`~pandas.DataFrame`
        A :class:`~pandas.DataFrame` object with model descriptions.
    """
    models_json = get_models_index()
    models = models_json["models"]
    filenames = [model['filename'] for model in models]
    descriptions = [model['details'] for model in models]
    return pd.DataFrame({'model': filenames, 'description': descriptions})
