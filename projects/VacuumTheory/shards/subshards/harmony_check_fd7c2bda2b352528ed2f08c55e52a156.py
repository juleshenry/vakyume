from cmath import *
from math import e, pi
import numpy as np


def check_harmony(A, Q_condensor_heat_duty, U, del_T_LM, **kwargs):
    return (A) - (Q_condensor_heat_duty / (U * del_T_LM))
