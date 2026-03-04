from cmath import *
from math import e, pi
import numpy as np


def check_harmony(Q_condensor_heat_duty, U_v, V_c, del_T_LM, **kwargs):
    return (V_c) - (Q_condensor_heat_duty / (U_v * del_T_LM))
