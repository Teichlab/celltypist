#CellTypist 1.0
from . import classifier, models, samples
from .annotate import annotate, Model
from .train import train
from .plot import dotplot
#CellTypist 2.0
from .contro.harmonize import harmonize
from .contro.integrate import integrate
from .contro.plot import treeplot
from .contro.plot import sankeyplot
from .contro.align import DistanceAlignment

__version__ = "1.4.0"
