from cmath import *
from math import e, pi
import numpy as np


def check_harmony(L_c, Q_condensor_heat_duty, c_p, del_T, **kwargs):
    return (L_c) - (Q_condensor_heat_duty / (c_p * del_T))
