import os
import numpy as np
from matplotlib import pyplot as plt

class Diamond_TL_Coherent():
    def __init__(self, ssp_file, bty_file, ati_file):
        self.ssp_file = ssp_file
        self.bty_file = bty_file
        self.ati_file = ati_file