import os
import sklearn
import celltypist
import unittest

# hide warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

BASE_PATH = os.path.dirname(__file__)
MODELS_PATH = os.path.join(BASE_PATH, "data", "models")
SAMPLES_PATH = os.path.join(BASE_PATH, "data", "samples")
INCLUDED_MODELS = ["brain", "default", "immune"]


class TestModels(unittest.TestCase):
    def test_list_models(self):
        self.assertTrue(INCLUDED_MODELS, celltypist.models.get_all())

    def test_load_defult_model(self):
        default_model = celltypist.models.load()
        self.assertIsInstance(default_model, sklearn.linear_model.SGDClassifier)

    def test_load_included_models(self):
        for model in INCLUDED_MODELS:
            named_model = celltypist.models.load(model)
            self.assertIsInstance(named_model, sklearn.linear_model.SGDClassifier)

    def test_load_custom_model(self):
        path_to_model = os.path.join(MODELS_PATH, "default.pkl")
        custom_model = celltypist.models.load(path_to_model)
        self.assertIsInstance(custom_model, sklearn.linear_model.SGDClassifier)

    def test_load_model_error(self):
        self.assertRaises(FileNotFoundError, celltypist.models.load, "asdasdasdasdasd")


class TestSamples(unittest.TestCase):
    def test_csv_filename(self):
        self.assertEqual(os.path.basename(celltypist.samples.csv()), "sample_cell_by_gene.csv")

    def test_csv_sample_content(self):
        with open(os.path.join(SAMPLES_PATH, "sample_cell_by_gene.csv")) as sample_file:
            test_sample = sample_file.read()
        with open(celltypist.samples.csv()) as package_file:
            package_sample = package_file.read()
        self.assertMultiLineEqual(test_sample, package_sample)
