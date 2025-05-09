import os
import unittest

import pandas as pd
from pathlib import Path
import pickle
import numpy as np
import warnings

warnings.filterwarnings("ignore")

from kmdcm.pydcm.dcm import (
    mdcm_set_up,
    eval_kernel,
)
from kmdcm.pydcm.dcm import FFE_PATH, espform, densform
from kmdcm.utils import dcm_utils as du
from kmdcm.pydcm.mdcm_dict import MDCM
from kmdcm.pydcm.kernel import KernelFit
cubes = list(Path("/pchem-data/meuwly/boittier/home/ref/comformation-cube/Benz").glob("esp*cube"))
str_cubes = [str(_) for _ in cubes]
model_dir = Path("/pchem-data/meuwly/boittier/home/mdcm_fast/mdcm/Benz/3-24")
mdcm_xyz = model_dir / "24_charges_refined.xyz"
mdcm_model = model_dir / "model.mdcm"
mdcm = mdcm_set_up(str_cubes, str_cubes, mdcm_cxyz=mdcm_xyz, mdcm_clcl=mdcm_model, local_pos=None)
