from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException
import numpy as np


class subshards:
    def check_harmony(H_1, H_2, P, P_P, **kwargs):
        return (P_P - P) - (H_2 - H_1)
