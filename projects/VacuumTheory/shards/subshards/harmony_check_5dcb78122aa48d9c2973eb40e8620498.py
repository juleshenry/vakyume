from cmath import *
from math import e, pi
import numpy as np


def check_harmony(NC, NS, SCON, installation_cost, **kwargs):
    return (installation_cost) - (16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35)
