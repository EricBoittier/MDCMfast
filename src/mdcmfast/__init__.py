from . import mdcm_fast, kmdcm
from .mdcm_fast import *
from .kmdcm import *
from pathlib import Path
import os

__all__ = ["mdcm_fast", "kmdcm"]

INSTALL_DIR = Path(os.path.abspath(__file__)).parent.parent.parent
print("Loaded MDCMFAST from", INSTALL_DIR)


def determine_n_charges(min_charges, max_charges):
    """
    Determine the number of charges to fit for a given molecule.
    """
    return range(min_charges, max_charges + 1)

