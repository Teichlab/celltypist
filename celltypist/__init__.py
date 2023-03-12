#CellTypist 1.0
from . import classifier, models, samples
from .annotate import annotate, Model
from .train import train
from .plot import dotplot
#CellTypist 2.0
from .contro.harmonize import harmonize, harmonise, DistanceAlignment
from .contro.integrate import integrate
from .contro.plot import treeplot, sankeyplot

__version__ = "1.4.0"
