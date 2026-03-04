from cmath import *
from math import e, pi
import numpy as np


def check_harmony(A, Q_condensor_heat_duty, U, del_T, **kwargs):
    return (Q_condensor_heat_duty) - (U * A * del_T)
